"""
PCA on highly variable genes
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from utils.processing import sc, USE_GPU


input_file = snakemake.input[0]
output_file = snakemake.output[0]
scale = snakemake.params.get('scale', False)
args = snakemake.params.get('args', {})
dask = snakemake.params.get('dask', False) and not USE_GPU and not scale
backed = snakemake.params.get('backed', False) and dask and not USE_GPU

logging.info(f'Read "{input_file}"...')
adata = read_anndata(
    input_file,
    var='var',
    obs='obs',
    uns='uns',
)
logging.info(adata.__str__())

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['scaled'] = scale
adata.uns['preprocessing']['pca'] = args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.obsm['X_pca'] = np.zeros((0, 30))
    adata.write_zarr(output_file)
    exit(0)

adata.X = read_anndata(
    input_file,
    X='X',
    dask=True,
    backed=True,
).X

# scaling TODO: move to separate rule
if scale:
    logging.info('Scale data...')
    adata = dask_compute(adata)
    sc.pp.scale(adata, max_value=10)

logging.info('Subset to highly variable genes...')
hvg_key = args.pop('mask_var', 'highly_variable')
adata, _ = subset_hvg(adata, var_column=hvg_key)

if USE_GPU:
    sc.get.anndata_to_GPU(adata)

logging.info(f'PCA with parameters: {args}')
sc.pp.pca(adata, **args)

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns', 'varm']
)
