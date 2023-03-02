"""
Build kNN graph on embedding
"""
import anndata
import scanpy as sc

input_file = snakemake.input.zarr
output_file = snakemake.output.zarr

adata = anndata.read_zarr(input_file)

sc.pp.neighbors(adata, use_rep='X_pca')

adata.write_zarr(output_file)
