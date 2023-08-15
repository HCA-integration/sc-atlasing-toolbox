from pathlib import Path
import shutil
import anndata as ad


def read_anndata(file):
    if file.endswith('.zarr'):
        adata = ad.read_zarr(file)
    else:
        adata = ad.read(file)
    return adata


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


def link_zarr(in_dir, out_dir, file_names=None, overwrite=False):
    """
    Link to existing zarr file
    """
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    
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
            print(f'Directory {new_file} exists, overwrite...')
            if new_file.is_dir():
                shutil.rmtree(new_file)
            else:
                new_file.unlink()
        new_file.symlink_to(f.resolve())