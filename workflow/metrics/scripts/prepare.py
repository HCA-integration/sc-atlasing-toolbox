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
from utils.processing import compute_neighbors, _filter_genes


def compute_pca(adata, matrix):
    X_pca, _, variance_ratio, variance = sc.tl.pca(matrix, return_info=True)
    adata.obsm['X_pca'] = X_pca
    adata.uns['pca'] = {
        'variance_ratio': variance_ratio,
        'variance': variance,
    }


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
    X=unintegrated_layer,
    obs='obs',
    obsm='obsm',
    obsp='obsp',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
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

# set HVGs
var_key = 'highly_variable' # TODO make configurable
new_var_column = 'metrics_features'
if var_key not in adata.var.columns:
    logging.info(f'{var_key} key not in adata var, setting all to True')
    adata.var[new_var_column] = True
else:
    adata.var[new_var_column] = adata.var[var_key]

logging.info('Filter all zero genes...')
all_zero_genes = _filter_genes(adata, min_cells=1)
adata.var[new_var_column] = adata.var[new_var_column] & ~adata.var_names.isin(all_zero_genes)

if output_type == 'full':
    hvg_matrix = subset_hvg(
        adata.copy(),
        var_column=new_var_column,
        compute_dask=True,
    )[0].X
    compute_pca(adata, matrix=hvg_matrix)
    del hvg_matrix
    force_neighbors = True
elif output_type == 'embed':
    logging.info('Run PCA on embedding...')
    compute_pca(adata, matrix=adata.obsm['X_emb'])

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
    files_to_keep=['obsm', 'obsp', 'raw', 'uns', 'var'],
    slot_map={'X': corrected_layer},
)

# unintegrated for comparison metrics
# if output_type != 'knn':  # assuming that there aren't any knn-based metrics that require PCA
logging.info(f'Prepare unintegrated data from layer={unintegrated_layer}...')
adata_raw = read_anndata(
    input_file,
    X=unintegrated_layer,
    dask=True,
    backed=True,
)

logging.info('Run PCA on unintegrated data...')
adata_raw.var = adata.var.copy()
hvg_matrix = subset_hvg(
    adata_raw.copy(),
    var_column=new_var_column,
    compute_dask=True,
)[0].X
compute_pca(adata_raw, matrix=hvg_matrix)
del hvg_matrix

raw_file = Path(output_file) / 'raw'
logging.info(f'Write to {raw_file}...')
write_zarr_linked(
    adata_raw,
    in_dir=output_file,
    out_dir=raw_file,
    files_to_keep=['layers', 'obsm', 'obsp', 'uns', 'var', 'varm', 'varp'],
    slot_map={'X': unintegrated_layer},
)