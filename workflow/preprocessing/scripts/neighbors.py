"""
Build kNN graph on embedding
"""
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

print('read...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

try:
    sc.pp.neighbors(adata, use_rep='X_pca', method='rapids')
except:
    print('Rapids failed, defaulting to UMAP implementation')
    sc.pp.neighbors(adata, use_rep='X_pca')

print('write...')
adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
