"""
PCA on highly variable genes
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError as e:
    logging.info('Importing rapids failed, using scanpy...')
    import scanpy as sc
    rapids = False
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, write_zarr_linked
from utils.misc import ensure_sparse


input_file = snakemake.input[0]
input_counts = snakemake.input.counts
output_file = snakemake.output[0]
scale = snakemake.params.get('scale', False)
args = snakemake.params.get('args', {})
dask = snakemake.params.get('dask', False) and not rapids and not scale
backed = snakemake.params.get('backed', False) and dask and not rapids

logging.info(f'Read "{input_file}"...')
adata = read_anndata(
    input_file,
    var='var',
    obs='obs',
    uns='uns',
)
logging.info(adata.__str__())

if adata.uns['preprocessing']['highly_variable_genes'] == False:
    args['use_highly_variable'] = False
else:
    args['use_highly_variable'] = True

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
    input_counts,
    X='X',
    var='var',
    backed=backed,
    dask=dask,
)[:, adata.var_names].X

if isinstance(adata.X, da.Array):
    logging.info('Convert dask chunks to sparse chunks...')
    adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype).persist()

# make sure data is on GPU for rapids_singlecell
if rapids:
    logging.info('Transfer to GPU...')
    sc.get.anndata_to_GPU(adata)

# scaling TODO: move to separate rule
if scale:
    logging.info('Scale data...')
    sc.pp.scale(adata, max_value=10)

logging.info('PCA...')
logging.info(f'parameters: {args}')
# args['n_comps'] = np.min([adata.n_obs-1, adata.n_vars-1, args.get('n_comps', 30)])
sc.pp.pca(adata, **args)

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns', 'varm']
)
