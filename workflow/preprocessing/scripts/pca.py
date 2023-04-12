"""
PCA on highly variable genes
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

sc.pp.pca(adata, use_highly_variable=True)
adata.obsm['X_pca'] = sparse.csr_matrix(adata.obsm['X_pca'])

print('write...')
adata.write(output_file, compression='lzf')
