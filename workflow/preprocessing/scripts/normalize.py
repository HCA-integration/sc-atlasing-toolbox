"""
Normalisation
"""
import anndata
import scanpy as sc

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

adata = anndata.read_zarr(input_zarr)
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.layers['normcounts'] = adata.X

adata.write_zarr(output_zarr)
