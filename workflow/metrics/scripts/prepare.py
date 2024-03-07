from pathlib import Path
import logging
logger = logging.getLogger('Preprocess for metrics')
logger.setLevel(logging.INFO)
import numpy as np
import anndata as ad

from utils.assertions import assert_pca
from utils.io import read_anndata, write_zarr_linked
from utils.processing import compute_neighbors, sc


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
label_key = params.get('label_key')
neighbor_args = params.get('neighbor_args', {})
unintegrated_layer = params.get('unintegrated_layer', 'X')
corrected_layer = params.get('corrected_layer', 'X')

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
    kwargs|= {'X': corrected_layer}
adata = read_anndata(input_file, **kwargs)

# remove cells without labels
n_obs = adata.n_obs
duplicated_mask = ~adata.obs.duplicated(keep='first')
adata = adata[duplicated_mask]  # remove duplicates
logger.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# logger.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logger.info(f'After: {adata.n_obs} cells')
if adata.is_view:
    adata = adata.copy()
force_neighbors = n_obs > adata.n_obs

logging.info(f'Computing neighbors for output type {output_type} force_neighbors={force_neighbors}...')
compute_neighbors(
    adata,
    output_type,
    force=force_neighbors,
    check_n_neighbors=True,
    **neighbor_args
)

logging.info(f'Prepare unintegrated data from layer={unintegrated_layer}...')
adata.raw = read_anndata(
    input_file,
    X=unintegrated_layer,
    var='var',
)[duplicated_mask]
logging.info('Run PCA on unintegrated data...')
adata_raw = adata.raw.to_adata()
sc.pp.pca(
    adata_raw,
    use_highly_variable='highly_variable' in adata_raw.var.columns
)
adata.obsm['X_pca'] = adata_raw.obsm['X_pca']
adata.varm['PCs'] = adata_raw.varm['PCs']
adata.uns['pca'] = adata_raw.uns['pca']
del adata_raw

# write to file
logging.info(adata.__str__())
logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'raw', 'var', 'varm', 'uns'],
    slot_map={
        'X': corrected_layer,
        'raw/X': unintegrated_layer,
    }
)
