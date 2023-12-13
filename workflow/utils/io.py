import warnings
import os
from pathlib import Path
import shutil
import anndata as ad
import zarr
import h5py
from scipy.sparse import csr_matrix
from anndata.experimental import read_elem, sparse_dataset


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


def read_anndata(
    file: str,
    dask: bool = False,
    backed: bool = False,
    **kwargs
) -> ad.AnnData:
    """
    Read anndata file
    :param file: path to anndata file in zarr or h5ad format
    :param kwargs: AnnData parameter to zarr group mapping
    """
    assert Path(file).exists(), f'File not found: {file}'
    
    if dask:
        read_func = read_dask
    elif backed:
        read_func = read_partial
    else:
        read_func = read_partial
    
    func, file_type = get_file_reader(file)
    
    f = func(file, 'r')
    kwargs = {x: x for x in f} if not kwargs else kwargs
    if len(f.keys()) == 0:
        return ad.AnnData()
    # check if keys are available
    for name, slot in kwargs.items():
        if slot not in f:
            warnings.warn(
                f'Cannot find "{slot}" for AnnData parameter `{name}` from "{file}"'
            )
    adata = read_func(f, backed=backed, **kwargs)
    if not backed and file_type == 'h5py':
        f.close()
    
    return adata


def read_partial(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    force_sparse_types: [str, list] = None,
    **kwargs
) -> ad.AnnData:
    """
    Partially read zarr or h5py groups
    :params group: file group
    :params force_sparse_types: encoding types to convert to sparse_dataset via csr_matrix
    :params backed: read sparse matrix as sparse_dataset
    :params **kwargs: dict of slot_name: slot, by default use all available slot for the zarr file
    :return: AnnData object
    """
    if force_sparse_types is None:
        force_sparse_types = []
    elif isinstance(force_sparse_types, str):
        force_sparse_types = [force_sparse_types]
    slots = {}
    for slot_name, slot in kwargs.items():
        print(f'Read slot "{slot}", store as "{slot_name}"...')
        if slot not in group:
            warnings.warn(f'Slot "{slot}" not found, skip...')
            slots[slot_name] = None
        else:
            elem = group[slot]
            iospec = ad._io.specs.get_spec(elem)
            if iospec.encoding_type in ("csr_matrix", "csc_matrix") and backed:
                slots[slot_name] = sparse_dataset(elem)
            elif iospec.encoding_type in force_sparse_types:
                slots[slot_name] = csr_matrix(read_elem(elem))
                if backed:
                    slots[slot_name] = sparse_dataset(slots[slot_name])
            else:
                slots[slot_name] = read_elem(elem)
    return ad.AnnData(**slots)


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
    from utils.sparse_dask import sparse_dataset_as_dask
    
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
            print('skip reading', elem_name)
            return None
        else:
            print('read', elem_name)
        
        if elem_name != '/' and all(elem_matches):
            print('skip reading', elem_name)
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
        print('Read as mudata...')
        return mu.read(file)
    elif file.endswith('.h5mu.zarr'):
        import mudata as mu
        print('Read as mudata from zarr...')
        return mu.read_zarr(file)
    else:
        print('Read as anndata...')
        return read_anndata(file)


def write_zarr(adata, file):
    def sparse_coo_to_csr(matrix):
        from dask.array import Array as DaskArray
        import sparse
        
        if isinstance(matrix, DaskArray) and isinstance(matrix._meta, sparse.COO):
            matrix = matrix.map_blocks(csr_matrix, dtype='float32')
        return matrix
    
    adata.X = sparse_coo_to_csr(adata.X)
    for layer in adata.layers:
        adata.layers[layer] = sparse_coo_to_csr(adata.layers[layer])
    adata.write_zarr(file) # doesn't seem to work with dask array


def link_zarr(in_dir, out_dir, file_names=None, overwrite=False, relative_path=True):
    """
    Link to existing zarr file
    """
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    
    if not in_dir.exists():
        return
    
    if file_names is None:
        file_names = [f.name for f in in_dir.iterdir()]

    if not out_dir.exists():
        out_dir.mkdir()
    for f in in_dir.iterdir():
        if f.name == '.snakemake_timestamp':
            continue  # skip snakemake timestamp
        if f.name not in file_names:
            continue
        new_file = out_dir / f.name
        if overwrite and new_file.exists():
            print(f'Replace {new_file} with link')
            if new_file.is_dir() and not new_file.is_symlink():
                shutil.rmtree(new_file)
            else:
                new_file.unlink()
        
        path_to_link_to = f.resolve()
        if relative_path:
            path_to_link_to = Path(
                os.path.relpath(
                    path_to_link_to,
                    new_file.parent.resolve()
                )
            )
        new_file.symlink_to(path_to_link_to)


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