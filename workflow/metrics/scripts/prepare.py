from pathlib import Path
import logging
logger = logging.getLogger('Preprocess for metrics')
logger.setLevel(logging.INFO)
import numpy as np
import anndata as ad
import scanpy as sc

from utils.assertions import assert_pca
from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked
from utils.processing import compute_neighbors


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
label_key = params.get('label_key')
neighbor_args = params.get('neighbor_args', {})
unintegrated_layer = params.get('unintegrated_layer', 'X')
corrected_layer = params.get('corrected_layer', 'X')

files_to_keep = ['obsm', 'obsp', 'var', 'uns']

# determine output types
# default output type is 'full'
output_type = read_anndata(input_file, uns='uns').uns.get('output_type', 'full')

logger.info(f'Read {input_file} ...')
kwargs = dict(
    obs='obs',
    obsm='obsm',
    obsp='obsp',
    var='var',
    varm='varm',
    uns='uns',
)
if output_type == 'full':
    kwargs |= {'X': corrected_layer}
adata = read_anndata(input_file, **kwargs)

# remove cells without labels
n_obs = adata.n_obs
# logger.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# logger.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logger.info(f'After: {adata.n_obs} cells')
# if adata.is_view:
#     adata = adata.copy()
force_neighbors = n_obs > adata.n_obs

if 'highly_variable' not in adata.var.columns:
    adata.var['highly_variable'] = True
adata, subsetted = subset_hvg(adata, var_column='highly_variable')
if subsetted:
    files_to_keep.extend(['X', 'varm', 'varp', 'raw'])

logging.info('Run PCA on integrated data...')
if output_type == 'full':
    sc.tl.pca(adata, mask_var='highly_variable')
elif output_type == 'embed':
    X_pca, _, variance_ratio, variance = sc.tl.pca(
        adata.obsm['X_emb'],
        return_info=True,
    )
    adata.obsm['X_pca'] = X_pca
    adata.uns['pca'] = {
        'variance_ratio': variance_ratio,
        'variance': variance,
    }
    print('after PCA', adata)

logging.info(f'Computing neighbors for output type {output_type} force_neighbors={force_neighbors}...')
compute_neighbors(
    adata,
    output_type,
    force=force_neighbors,
    check_n_neighbors=True,
    **neighbor_args
)

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
    slot_map={'X': corrected_layer},
)

# unintegrated for comparison metrics
if output_type != 'knn':  # assuming that there aren't any knn-based metrics that require PCA
    logging.info(f'Prepare unintegrated data from layer={unintegrated_layer}...')
    adata_raw = read_anndata(
        input_file,
        X=unintegrated_layer,
        var='var',
        dask=True,
        backed=True,
    )

    if 'highly_variable' not in adata_raw.var.columns:
        adata_raw.var['highly_variable'] = True
    adata_raw, subsetted = subset_hvg(adata_raw, var_column='highly_variable')
    if subsetted:
        files_to_keep.extend(['X', 'varm', 'varp'])

    logging.info('Run PCA on unintegrated data...')
    sc.tl.pca(adata_raw, mask_var='highly_variable')

    logging.info(f'Write to {output_file}/raw...')
    write_zarr_linked(
        adata_raw,
        in_dir=input_file,
        out_dir=f'{output_file}/raw',
        files_to_keep=files_to_keep,
        slot_map={'X': unintegrated_layer},
    )