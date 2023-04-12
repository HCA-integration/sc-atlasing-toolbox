"""
UMAP
"""
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

print('read...')
adata = read_anndata(input_file)

try:
    sc.tl.umap(adata, method='rapids')
except:
    print('Rapids failed, defaulting to UMAP implementation')
    sc.tl.umap(adata)
adata.obsm['X_umap'] = sparse.csr_matrix(adata.obsm['X_umap'])

print('write...')
adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
