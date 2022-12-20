"""
PCA on highly variable genes
"""
import anndata
import scanpy as sc

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

adata = anndata.read_zarr(input_zarr)

sc.pp.pca(adata)

adata.write_zarr(output_zarr)
