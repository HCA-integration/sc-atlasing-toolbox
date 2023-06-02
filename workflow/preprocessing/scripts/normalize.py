"""
Normalisation
"""
from scipy import sparse
import scanpy as sc

from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

print('read...')
adata = read_anndata(input_file)
adata.X = sparse.csr_matrix(adata.X)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.X = sparse.csr_matrix(adata.X)

adata.layers['normcounts'] = adata.X

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = 'default'
adata.uns['preprocessing']['log-transformed'] = True

print('write...')
adata.write(output_file, compression='lzf')
