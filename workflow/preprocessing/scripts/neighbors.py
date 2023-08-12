"""
Build kNN graph on embedding
"""
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
from utils.io import read_anndata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = {k: v for k, v in snakemake.params.args.items()}

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

# set defaults
if 'use_rep' not in params:
    if 'X_pca' in adata.obsm:
        params['use_rep'] = 'X_pca'
        params['n_pcs'] = adata.obsm['X_pca'].shape[1]
    else:
        params['use_rep'] = 'X'
logging.info(str(params))

# compute kNN graph
try:
    logging.info('Compute kNN graph...')
    sc.pp.neighbors(adata, method='rapids', **params)
except Exception as e:
    logging.info(e)
    logging.info('Rapids failed, defaulting to UMAP implementation')
    sc.pp.neighbors(adata, **params)

# remove redundant data
del adata.obsm['X_pca']

logging.info(f'Write to {output_file}...')
# if not adata.uns['preprocessing']['scaled']:
#     adata.X = sparse.csr_matrix(adata.X)
adata.write_zarr(output_file)
