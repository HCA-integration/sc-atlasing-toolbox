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

from utils.io import read_anndata, write_zarr_linked, link_zarr
from utils.misc import ensure_sparse


input_dir = snakemake.input[0]
output_dir = snakemake.output[0]

logging.info(f'Read {input_dir}...')
adata = read_anndata(
    input_dir,
    X=snakemake.params.get('raw_counts', 'X'),
    obs='obs',
    var='var',
    uns='uns'
)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.X = np.zeros(adata.shape)
    adata.write_zarr(output_dir)
    exit(0)

if input_dir.endswith('.h5ad'):
    adata.layers['counts'] = adata.X

# make sure data is on GPU for rapids_singlecell
if rapids:
    adata.X = adata.X.astype('float32')
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

logging.info(f'Write to {output_dir}...')
if input_dir.endswith('.h5ad'):
    adata.layers['normcounts'] = adata.X

logging.info(f'Write to {output_dir}...')
write_zarr_linked(
    adata,
    input_dir,
    output_dir,
    files_to_keep=['X', 'uns'],
    slot_map={
        'raw': '',
        'layers/counts': 'X',
    }
)

if input_dir.endswith('.zarr'):
    logging.info('Linking layers...')
    output_dir = Path(output_dir)
    link_zarr(
        in_dir=output_dir / 'X',
        out_dir=output_dir / 'layers' / 'normcounts',
        overwrite=True
    )