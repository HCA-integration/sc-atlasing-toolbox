from pathlib import Path
import logging
logger = logging.getLogger('Preprocess for metrics')
logger.setLevel(logging.INFO)
import scanpy as sc

from utils.assertions import assert_pca
from utils.accessors import adata_to_memory
from utils.io import read_anndata, write_zarr_linked
from utils.processing import compute_neighbors


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
label_key = params.get('label_key')
neighbor_args = params.get('neighbor_args', {})
unintegrated_layer = params.get('unintegrated_layer')

files_to_keep = ['obsp', 'var', 'uns']

logging.info('Read file...')
adata = read_anndata(
    input_file,
    X=unintegrated_layer,
    obs='obs',
    var='var',
    obsm='obsm',
    obsp='obsp',
    varm='varm',
    uns='uns',
    backed=True,
)

# determine output types
# default output type is 'full'
output_type = adata.uns.get('output_type', 'full')

if output_type == 'full':
    adata.layers['corrected'] = read_anndata(input_file, X='X', backed=True).X

try:
    assert_pca(adata)
except AssertionError as e:
    logging.info(f'PCA not found: {e}\nrecomputing PCA...')
    adata = adata_to_memory(adata)
    sc.pp.pca(adata, mask='highly_variable')
    files_to_keep.extend(['obsm', 'varm'])


n_obs = adata.n_obs
# logger.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# TODO: Error when connectivities/distances have different dimensions than X
# logger.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logger.info(f'After: {adata.n_obs} cells')
force_neighbors = n_obs > adata.n_obs

logging.info(f'Computing neighbors for output type {output_type}...')
compute_neighbors(adata, output_type, force=force_neighbors, **neighbor_args)

# TODO: compute for unintegrated
# if output_type == 'full':
#     adata.X = adata.layers['corrected']
# neighbor_args['use_rep'] = 'X'
# compute_neighbors(adata, 'full', force=force_neighbors, **neighbor_args)

# write to file
logging.info(adata.__str__())
logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
)
