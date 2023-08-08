"""
Build kNN graph on embedding
"""
import logging
logging.basicConfig(level=logging.INFO)
import sys
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.args

if args is None:
    args = {}
logging.info(str(args))

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

try:
    logging.info('Compute kNN graph...')
    sc.pp.neighbors(adata, method='rapids', **args)
except Exception as e:
    logging.info(e)
    logging.info('Rapids failed, defaulting to UMAP implementation')
    sc.pp.neighbors(adata, **args)

# remove redundant data
del adata.obsm['X_pca']

logging.info(f'Write to {output_file}...')
# if not adata.uns['preprocessing']['scaled']:
#     adata.X = sparse.csr_matrix(adata.X)
adata.write_zarr(output_file)
