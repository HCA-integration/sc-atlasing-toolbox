"""
PCA on highly variable genes
"""
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
scale = snakemake.params['scale']

print('read...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

# scaling
if scale:
    sc.pp.scale(adata, max_value=10, zero_center=True)
    sc.pp.pca(adata, use_highly_variable=True)
else:
    sc.pp.pca(adata, use_highly_variable=True)

adata.obsm['X_pca'] = sparse.csr_matrix(adata.obsm['X_pca'])

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['scaled'] = scale

print('write...')
adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
