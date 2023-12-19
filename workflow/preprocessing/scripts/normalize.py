"""
Normalisation
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    rapids = False

from utils.io import read_anndata, link_zarr_partial
from utils.misc import ensure_sparse


input_file = snakemake.input[0]
output_file = snakemake.output[0]

# select counts layer
logging.info('Select layer...')
layer = snakemake.params.get('raw_counts', 'X')
if layer is None:
    layer = 'X'

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, X=layer, obs='obs', var='var', uns='uns')

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.X = np.zeros(adata.shape)
    adata.write_zarr(output_file)
    exit(0)

# make sure data is on GPU for rapids_singlecell
if rapids:
    sc.utils.anndata_to_GPU(adata)

logging.info('normalize_total...')
sc.pp.normalize_total(adata)
logging.info('log-transform...')
sc.pp.log1p(adata)

if rapids:
    sc.utils.anndata_to_CPU(adata)

ensure_sparse(adata)

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = 'default'
adata.uns['preprocessing']['log-transformed'] = True

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['X', 'uns'])