import warnings
import os
from pathlib import Path
import shutil
import anndata as ad
import zarr


def read_anndata(file, **kwargs):
    """
    Read anndata file
    :param file: path to anndata file
    :param kwargs: kwargs for partial zarr reader
    """
    if file.endswith('.zarr'):
        adata = read_zarr_partial(file, **kwargs)
    elif file.endswith('.h5ad'):
        adata = ad.read_h5ad(file)
    else:
        raise ValueError(f'Unknown file format: {file}')
    return adata


def read_zarr_partial(file, **kwargs):
    """
    Partially read zarr files
    :param file: path to zarr file
    :param kwargs: dict of slot_name: slot, by default use all available slot for the zarr file
    :return: AnnData object
    """
    slots = {}
    with zarr.open(file) as z:
        if not kwargs:
            kwargs = {x: x for x in z}
        for slot_name, slot in kwargs.items():
            print(f'Read slot "{slot}", store as "{slot_name}"...')
            if slot not in z:
                warnings.warn(f'Slot "{slot}" not found in zarr file, skip...')
                slots[slot_name] = None
            else:
                slots[slot_name] = ad.experimental.read_elem(z[slot])
        return ad.AnnData(**slots)


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