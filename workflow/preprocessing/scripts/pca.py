"""
PCA on highly variable genes
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    rapids = False

from utils.io import read_anndata, link_zarr
from utils.misc import ensure_sparse


input_file = snakemake.input[0]
input_counts = snakemake.input.counts
output_file = snakemake.output[0]
scale = snakemake.params.get('scale', False)
args = snakemake.params.get('args', {})

logging.info(f'Read "{input_file}"...')
adata = read_anndata(input_file, var='var', obs='obs', uns='uns')

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

adata.X = read_anndata(input_counts, X='X', var='var', obs='obs')[:, adata.var_names].X
ensure_sparse(adata)

# make sure data is on GPU for rapids_singlecell
if rapids:
    sc.utils.anndata_to_GPU(adata)

# scaling
if scale:
    logging.info('Scale data...')
    sc.pp.scale(adata, max_value=10)

logging.info('PCA...')
# args['n_comps'] = np.min([adata.n_obs-1, adata.n_vars-1, args.get('n_comps', 30)])
sc.pp.pca(adata, use_highly_variable=True, n_comps=50, **args)

logging.info(f'Write to "{output_file}"...')
del adata.X
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_link = [f for f in input_files if f not in ['obsm', 'uns', 'varm']]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_link,
        overwrite=True,
    )
