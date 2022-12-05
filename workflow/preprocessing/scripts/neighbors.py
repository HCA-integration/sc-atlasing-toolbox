"""
Build kNN graph on embedding
"""
import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad

adata = sc.read(input_h5ad)

sc.pp.neighbors(adata, use_rep='X_pca')

adata.write(output_h5ad)
