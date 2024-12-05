import warnings
import os
from pathlib import Path
from typing import MutableMapping
import shutil
import json
import anndata as ad
import zarr
import h5py
import numpy as np
from scipy.sparse import csr_matrix
from anndata.experimental import read_elem, sparse_dataset
from functools import partial
from dask import array as da


def print_flushed(*args, **kwargs):
    kwargs |= dict(flush=True)
    verbose = kwargs.pop('verbose', True)
    if verbose:
        print(*args, **kwargs)


zarr.default_compressor = zarr.Blosc(shuffle=zarr.Blosc.SHUFFLE)


def get_file_reader(file):
    if file.endswith(('.zarr', '.zarr/')):
        func = zarr.open
        file_type = 'zarr'
    elif file.endswith('.h5ad'):
        func = h5py.File
        file_type = 'h5py'
    else:
        raise ValueError(f'Unknown file format: {file}')
    return func, file_type


def check_slot_exists(file, slot):
    func, file_type = get_file_reader(file)
    with func(file, 'r') as f:
        exists = slot in f
    return exists


def to_memory(matrix):
    if isinstance(matrix, (ad.experimental.CSRDataset, ad.experimental.CSCDataset)):
        print_flushed('Convert to memory...')
        return matrix.to_memory()
    return matrix


def read_anndata(
    file: str,
    dask: bool = False,
    backed: bool = False,
    fail_on_missing: bool = True,
    exclude_slots: list = None,
    chunks: [int, tuple] = ('auto', -1),
    stride: int = 1000,
    **kwargs
) -> ad.AnnData:
    """
    Read anndata file
    :param file: path to anndata file in zarr or h5ad format
    :param kwargs: AnnData parameter to zarr group mapping
    """
    # assert Path(file).exists(), f'File not found: {file}'
    if exclude_slots is None:
        exclude_slots = []
    elif exclude_slots == 'all':
        exclude_slots = ['X', 'layers', 'raw']

    func, file_type = get_file_reader(file)
    try:
        store = func(file, 'r')
    except zarr.errors.PathNotFoundError as e:
        raise FileNotFoundError(f'Cannot read file {file}') from e
    
    # set default kwargs
    kwargs = {x: x for x in store} if not kwargs else kwargs
    # set key == value if value is None
    kwargs |= {k: k for k, v in kwargs.items() if v is None}
    # exclude slots
    kwargs = {k: v for k, v in kwargs.items() if k not in exclude_slots}
    
    # return an empty AnnData object if no keys are available
    if len(store.keys()) == 0:
        return ad.AnnData()
    
    # check if keys are available
    for name, slot in kwargs.items():
        if slot not in store:
            message = f'Cannot find "{slot}" for AnnData parameter `{name}`'
            message += f'\nfile: {file}\navailable slots: {list(store)}'
            if fail_on_missing:
                raise ValueError(message)
            warnings.warn(f'{message}, will be skipped')
    adata = read_partial(
        store,
        dask=dask,
        backed=backed,
        chunks=chunks,
        stride=stride,
        **kwargs
    )
    if not backed and file_type == 'h5py':
        store.close()
    
    return adata


def read_partial(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    dask: bool = False,
    chunks: [int, tuple] = ('auto', -1),
    stride: int = 1000,
    force_sparse_types: [str, list] = None,
    force_sparse_slots: [str, list] = None,
    verbose: bool = False,
    **kwargs
) -> ad.AnnData:
    """
    Partially read zarr or h5py groups
    :params group: file group
    :params force_sparse_types: encoding types to convert to sparse_dataset via csr_matrix
    :params backed: read sparse matrix as sparse_dataset
    :params dask: read any matrix as dask array
    :params chunks: chunks parameter for creating dask array
    :params stride: stride parameter for creating backed dask array
    :params **kwargs: dict of to_slot: slot, by default use all available slot for the zarr file
    :return: AnnData object
    """
    if force_sparse_types is None:
        force_sparse_types = []
    elif isinstance(force_sparse_types, str):
        force_sparse_types = [force_sparse_types]
        
    if force_sparse_slots is None:
        force_sparse_slots = []
    elif isinstance(force_sparse_slots, str):
        force_sparse_slots = [force_sparse_slots]
    force_sparse_slots.extend(['X', 'layers/', 'raw/X'])
    
    print_flushed(f'dask: {dask}, backed: {backed}')
    if dask:
        print_flushed('chunks:', chunks)
    
    slots = {}
    for to_slot, from_slot in kwargs.items():
        print_flushed(f'Read slot "{from_slot}", store as "{to_slot}"...', verbose=verbose)
        force_slot_sparse = any(from_slot.startswith((x, f'/{x}')) for x in force_sparse_slots)
        if from_slot in ['layers', '/layers', 'raw', '/raw']:
            slots[to_slot] = {
                sub_slot: read_slot(
                    group,
                    f'{from_slot}/{sub_slot}',
                    force_sparse_types,
                    force_slot_sparse,
                    backed=backed,
                    dask=dask,
                    chunks=chunks,
                    stride=stride,
                    verbose=verbose,
                )
                for sub_slot in group[from_slot]
            }
        else:
            slots[to_slot] = read_slot(
                group,
                from_slot,
                force_sparse_types,
                force_slot_sparse,
                backed=backed,
                dask=dask,
                chunks=chunks,
                stride=stride,
                verbose=verbose,
            )

    return ad.AnnData(**slots)


def read_slot(
    group: [h5py.Group, zarr.Group],
    slot: str,
    force_sparse_types: list,
    force_slot_sparse: bool,
    backed: bool,
    dask: bool,
    chunks: [int, tuple],
    stride: int,
    verbose: bool,
):
    if slot not in group:
        warnings.warn(f'Slot "{slot}" not found, skip...')
        return None
    
    if dask:
        return _read_slot_dask(
            group,
            slot,
            force_sparse_types,
            force_slot_sparse,
            stride=stride,
            chunks=chunks,
            backed=backed,
            verbose=verbose,
        )
    return _read_slot_default(
        group,
        slot,
        force_sparse_types,
        force_slot_sparse,
        backed=backed,
        verbose=verbose,
    )


def _read_slot_dask(
    group,
    slot,
    force_sparse_types,
    force_slot_sparse,
    stride,
    chunks,
    backed,
    verbose,
):
    """
    :param group: h5py.Group or zarr.Group
    :param slot: slot name
    :stride: stride parameter for creating backed dask array, ignored when backed=False
    :chunks: chunks parameter for creating dask array ignored when backed=True
    """
    from .sparse_dask import sparse_dataset_as_dask, read_as_dask_array
    
    elem = group[slot]
    iospec = ad._io.specs.get_spec(elem)
    
    if iospec.encoding_type in ("csr_matrix", "csc_matrix"):
        if backed:
            print_flushed(f'Read {slot} as backed sparse dask array...', verbose=verbose)
            elem = sparse_dataset(elem)
            return sparse_dataset_as_dask(elem, stride=stride)
        print_flushed(f'Read {slot} as sparse dask array...', verbose=verbose)
        elem = read_as_dask_array(elem, chunks=chunks)
        return elem.map_blocks(csr_matrix_int64_indptr, dtype=elem.dtype)
    
    elif iospec.encoding_type in force_sparse_types or force_slot_sparse:
        print_flushed(f'Read {slot} as dask array and convert blocks to csr_matrix...', verbose=verbose)
        elem = read_as_dask_array(elem, chunks=chunks)
        return elem.map_blocks(csr_matrix_int64_indptr, dtype=elem.dtype)
    
    elif iospec.encoding_type == "array":
        print_flushed(f'Read {slot} as dask array...', verbose=verbose)
        return read_as_dask_array(elem, chunks=chunks)
    
    return read_elem(elem)


def _read_slot_default(group, slot, force_sparse_types, force_slot_sparse, backed, verbose):
    elem = group[slot]
    iospec = ad._io.specs.get_spec(elem)
    
    if iospec.encoding_type in ("csr_matrix", "csc_matrix"):
        if backed:
            print_flushed(f'Read {slot} as backed sparse matrix...', verbose=verbose)
            return sparse_dataset(elem)
        return read_elem(elem)
    elif iospec.encoding_type in force_sparse_types or force_slot_sparse:
        print_flushed(f'Read {slot} and convert to csr matrix...', verbose=verbose)
        return csr_matrix(read_elem(elem))
    else:
        return read_elem(elem)


def csr_matrix_int64_indptr(x):
    x = csr_matrix(x)
    x.indptr = x.indptr.astype(np.int64)
    x.indices = x.indices.astype(np.int64) # seems to be necessary to avoid "ValueError: Output dtype not compatible with inputs."
    return x


# deprecated
def read_dask(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    obs_chunk: int = 1000,
    **kwargs
) -> ad.AnnData:
    """
    Modified from https://anndata.readthedocs.io/en/latest/tutorials/notebooks/%7Bread%2Cwrite%7D_dispatched.html
    """
    from anndata.experimental import read_dispatched
    from .sparse_dask import sparse_dataset_as_dask
    
    def callback(func, elem_name: str, elem, iospec):
        import re
        import dask.array as da
        import sparse
        
        elem_matches = [
            not (
                bool(re.match(f'/{e}(/.?|$)', elem_name)) 
                or f'/{e}'.startswith(elem_name) 
            )
            for e in kwargs.values()
        ]
        if elem_name != '/' and all(elem_matches):
            print_flushed('skip reading', elem_name)
            return None
        else:
            print_flushed('read', elem_name)
        
        if elem_name != '/' and all(elem_matches):
            print_flushed('skip reading', elem_name)
            return None
        elif iospec.encoding_type in (
            "dataframe",
            "awkward-array",
        ):
            # Preventing recursing inside of these types
            return read_elem(elem)
        elif iospec.encoding_type in ("csr_matrix", "csc_matrix"):
            # return da.from_array(read_elem(elem))
            matrix = sparse_dataset_as_dask(sparse_dataset(elem), obs_chunk)
            return matrix.map_blocks(sparse.COO)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        return func(elem)

    return read_dispatched(group, callback=callback)


# deprecated
def read_anndata_or_mudata(file):
    if file.endswith('.h5mu'):
        import mudata as mu
        print_flushed('Read as mudata...')
        return mu.read(file)
    elif file.endswith('.h5mu.zarr'):
        import mudata as mu
        print_flushed('Read as mudata from zarr...')
        return mu.read_zarr(file)
    else:
        print_flushed('Read as anndata...')
        return read_anndata(file)


def write_zarr(adata, file):
    def sparse_coo_to_csr(matrix):
        from dask import array as da
        import sparse
        
        if isinstance(matrix, da.Array) and isinstance(matrix._meta, sparse.COO):
            matrix = matrix.map_blocks(lambda x: x.tocsr(), dtype='float32')
        return matrix
    
    adata.X = sparse_coo_to_csr(adata.X)
    for layer in adata.layers:
        adata.layers[layer] = sparse_coo_to_csr(adata.layers[layer])
    adata.write_zarr(file) # doesn't seem to work with dask array


def link_file(in_file, out_file, relative_path=True, overwrite=False):
    in_file = Path(in_file).resolve(True)
    out_file = Path(out_file)
    out_dir = out_file.parent.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    
    if relative_path:
        in_file = Path(os.path.relpath(in_file, out_dir))
    
    if overwrite and out_file.exists():
        if out_file.is_dir() and not out_file.is_symlink():
            print_flushed(f'replace {out_file}...')
            shutil.rmtree(out_file)
        else:
            out_file.unlink()
    
    out_file.symlink_to(in_file)
    assert out_file.exists(), f'Linking failed for {out_file.resolve()} -> {in_file}'


def link_zarr(
    in_dir: [str, Path],
    out_dir: [str, Path],
    file_names: list = None,
    overwrite: bool = False,
    relative_path: bool = True,
    slot_map: MutableMapping = None,
    in_dir_map: MutableMapping = None,
):
    """
    Link to existing zarr file
    :param in_dir: path to existing zarr file or mapping of input slot in slot_map to input path
    :param out_dir: path to output zarr file
    :param file_names: list of files to link, if None, link all files
    :param overwrite: overwrite existing output files
    :param relative_path: use relative path for link
    :param slot_map: custom mapping of output slot to input slots
    :param in_dir_map: input directory map for input slots
    :param kwargs: custom mapping of output slot to input slot,
        will update default mapping of same input and output naming
    """
    def prune_nested_links(slot_map, in_dir_map):
        # determine equivalence classes of slots (top hierarchy)
        eq_classes = {}
        for out_slot in slot_map:
            if '/' not in out_slot:
                continue
            eq = out_slot.rsplit('/', 1)[0]
            if eq not in slot_map:
                continue
            eq_classes.setdefault(eq, []).append(out_slot)

        for out_slot in eq_classes.keys():
            in_slot = slot_map[out_slot]
            # link all files of the equivalence class
            for f in (in_dir_map[in_slot] / in_slot).iterdir():
                new_out_slot = f'{out_slot}/{f.name}'
                new_in_slot = f'{in_slot}/{f.name}'
                # skip if already specified or .snakefile_timestamp
                if new_out_slot in slot_map or f.name == '.snakemake_timestamp':
                    continue
                # inherit in_dir
                in_dir_map[new_in_slot] = in_dir_map[in_slot]
                # update slot_map
                slot_map[new_out_slot] = new_in_slot
            # remove equivalence class once done to avoid overwriting
            del slot_map[out_slot]
        return slot_map, in_dir_map
    
    if file_names is None:
        file_names = [] if in_dir is None else [f.name for f in Path(in_dir).iterdir()]
    file_names = [
        file for file in file_names
        if file not in ('.snakemake_timestamp')
    ]
    
    if slot_map is None:
        slot_map = {}
    
    slot_map = {file: file for file in file_names} | slot_map
    slot_map |= {k: k for k, v in slot_map.items() if v is None}
    
    if in_dir_map is None:
        in_dir_map = {}
    
    in_dir_map = {slot: in_dir for slot in slot_map.values()} | in_dir_map
    in_dir_map = {slot: Path(path) for slot, path in in_dir_map.items()}
    for _dir in in_dir_map.values():
        assert _dir.exists(), f'Input directory {_dir} does not exist'
    
    # deal with nested mapping
    slot_map, in_dir_map = prune_nested_links(slot_map, in_dir_map)

    # link all files
    out_dir = Path(out_dir)
    slot_map = sorted(
        slot_map.items(),
        key=lambda item: out_dir.name in str(in_dir_map[item[1]]),
        reverse=False,
    )
    print_flushed('slot_map:', pformat(slot_map))

    for out_slot, in_slot in slot_map:
        in_dir = in_dir_map[in_slot]
        in_file_name = str(in_dir).split('.zarr')[-1] + '/' + in_slot
        out_file_name = str(out_dir).split('.zarr')[-1] + '/' + out_slot
        print_flushed(f'Link {out_file_name} -> {in_file_name}')
        link_file(
            in_file=in_dir / in_slot,
            out_file=out_dir / out_slot,
            relative_path=relative_path,
            overwrite=overwrite,
        )


# TODO: deprecate
def link_zarr_partial(in_dir, out_dir, files_to_keep=None, overwrite=True, relative_path=True):
    """
    Link zarr files excluding defined slots
    """
    if not in_dir.endswith('.zarr'):
        return
    if files_to_keep is None:
        files_to_keep = []
    in_dirs = [f.name for f in Path(in_dir).iterdir()]
    files_to_link = [f for f in in_dirs if f not in files_to_keep]
    link_zarr(
        in_dir=in_dir,
        out_dir=out_dir,
        file_names=files_to_link,
        overwrite=overwrite,
        relative_path=relative_path,
    )


def write_zarr_linked(
    adata: ad.AnnData,
    in_dir: [str, Path],
    out_dir: [str, Path],
    relative_path: bool = True,
    files_to_keep: list = None,
    slot_map: MutableMapping = None,
    in_dir_map: MutableMapping = None,
):
    """
    Write adata to linked zarr file
    :param adata: AnnData object
    :param in_dir: path to existing zarr file
    :param out_dir: path to output zarr file
    :param files_to_keep: list of files to keep and not overwrite
    :param relative_path: use relative path for link
    :param slot_map: custom mapping of output slot to input slot, for slots that are not in files_to_keep
    """
    if in_dir is None:
        in_dirs = []
    else:
        in_dir = Path(in_dir)
        if not in_dir.name.endswith('.zarr'):
            adata.write_zarr(out_dir)
            return
        in_dirs = [f.name for f in in_dir.iterdir()]
    
    if files_to_keep is None:
        files_to_keep = []
    
    # Unique the list
    files_to_keep = list(set(files_to_keep))

    # Make a list of existing subpaths to keep                    
    file_to_link_clean = [
        x.split('/')[0] if not x.startswith('/')  # take top level directory
        else x.split('/')[1]  # take first directory after root /
        for x in files_to_keep
    ]

    # For those not keeping, link
    files_to_link = [
        f for f in in_dirs
        if f.split('/', 1)[-1] not in file_to_link_clean
    ]
    
    if slot_map is None:
        slot_map = {}
    extra_slots_to_link = list(slot_map.keys())
    
    # keep only slots that are not explicitly in files_to_keep
    slot_map = {
        in_slot: out_slot 
        for in_slot, out_slot in slot_map.items()
        if in_slot not in files_to_keep
    }
    
    # remove slots that will be overwritten anyway
    for slot in set(files_to_link+extra_slots_to_link):
        if hasattr(adata, slot):
            print_flushed(f'remove slot to be linked: {slot}')
            delattr(adata, slot)
    
    # write zarr file
    adata.write_zarr(out_dir)
    
    # link files
    link_zarr(
        in_dir=in_dir,
        out_dir=out_dir,
        file_names=files_to_link,
        overwrite=True,
        relative_path=relative_path,
        slot_map=slot_map,
        in_dir_map=in_dir_map,
    )
    
    # update .zattrs files
    for slot in ['obs', 'var']:
        zattrs_file = Path(out_dir) / slot / '.zattrs'

        if not (zattrs_file).exists() or slot in file_to_link_clean:
            continue
        
        with open(zattrs_file, 'r') as file:
            zattrs = json.load(file)
        
        # add all columns (otherwise linked columns are not included)
        columns = [
            f.name for f in zattrs_file.parent.iterdir()
            if not f.name.startswith('.') and f.name != zattrs['_index']
        ]
        zattrs['column-order'] = sorted(columns)
        
        with open(zattrs_file, 'w') as file:
            json.dump(zattrs, file, indent=4)
