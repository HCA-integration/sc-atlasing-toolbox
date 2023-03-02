"""
UMAP
"""
import anndata
import scanpy as sc

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

adata = anndata.read_zarr(input_zarr)

sc.tl.umap(adata)

adata.write_zarr(output_zarr)
