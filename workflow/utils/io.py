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
        adata = mu.read(file)
    else:
        print('Read as anndata...')
        adata = read_anndata(file)
    return adata