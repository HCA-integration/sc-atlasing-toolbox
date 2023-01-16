import anndata
import scanpy as sc

# TODO: put in common location
def read_anndata(file):
    if file.endswith('.zarr'):
        adata = anndata.read_zarr(file)
    else:
        adata = sc.read(file)
    return adata
