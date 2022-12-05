"""
PCA on highly variable genes
"""
import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad

adata = sc.read(input_h5ad)

sc.pp.pca(adata)

adata.write(output_h5ad)
