"""
Build kNN graph on embedding
"""
import sys
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.args

if args is None:
    args = {}
print(args)

print('read...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

try:
    sc.pp.neighbors(adata, method='rapids', **args)
except:
    print('Rapids failed, defaulting to UMAP implementation', file=sys.stderr)
    sc.pp.neighbors(adata, **args)

print('write...')
if not adata.uns['preprocessing']['scaled']:
    adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
