"""
UMAP
"""
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

try:
    logging.info('Compute UMAP...')
    sc.tl.umap(adata, method='rapids')
except:
    logging.info('Rapids failed, defaulting to UMAP implementation')
    sc.tl.umap(adata)
# adata.obsm['X_umap'] = sparse.csr_matrix(adata.obsm['X_umap'])

logging.info(f'Write to {output_file}...')
if not adata.uns['preprocessing']['scaled']:
    adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
