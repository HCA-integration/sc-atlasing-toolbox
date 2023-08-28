"""
Normalisation
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
adata.X = sparse.csr_matrix(adata.X)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

# adata.layers['counts'] = adata.X.copy()
# select counts layer
logging.info('Select layer...')
layer = snakemake.params.get('raw_counts', 'X')
adata.X = adata.X if layer == 'X' or layer is None else adata.layers[layer]

logging.info('normalize_total...')
sc.pp.normalize_total(adata)
logging.info('log-transform...')
sc.pp.log1p(adata)
adata.X = sparse.csr_matrix(adata.X)

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = 'default'
adata.uns['preprocessing']['log-transformed'] = True

logging.info(f'Write to {output_file}...')
del adata.raw
del adata.layers
adata.write_zarr(output_file)
