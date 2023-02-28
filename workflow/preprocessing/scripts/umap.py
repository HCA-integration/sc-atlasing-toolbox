"""
UMAP
"""
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

print('read...')
adata = read_anndata(input_file)

sc.tl.umap(adata, method='rapids')

print('write...')
adata.write(output_file)
