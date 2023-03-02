"""
Build kNN graph on embedding
"""
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

print('read...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

sc.pp.neighbors(adata, use_rep='X_pca', method='rapids')

print('write...')
adata.write(output_file)
