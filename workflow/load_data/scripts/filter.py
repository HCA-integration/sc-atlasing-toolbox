"""
Data filtering
"""
from pprint import pformat
import numpy as np
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils_pipeline.io import read_anndata, link_zarr, to_memory


def save_empty(input_zarr, output_zarr, output_removed):
    # link unfiltered
    link_zarr(
        in_dir=input_zarr,
        out_dir=output_zarr,
    )

    # save removed cells
    adata_empty = read_anndata(input_zarr, var='var', uns='uns')
    adata_empty.X = np.zeros(shape=(0, adata_empty.var.shape[0]))
    adata_empty.write_zarr(output_removed)
    exit(0)


input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
output_removed = snakemake.output.removed
params = snakemake.params.get('filter', {})

logging.info(pformat(params))

if len(params) == 0:
    save_empty(input_zarr, output_zarr, output_removed)

backed = params.get('backed', True)
adata = read_anndata(input_zarr, backed=backed)
logging.info(adata.__str__())

keep_indices = []

# subset to available annotations TODO: make this more generic
if 'cell_type' in adata.obs.columns:
    keep_indices.extend(
        adata[adata.obs['cell_type'].notna()].obs_names
    )

# explicit filter
if 'remove_by_colum' in params:
    logging.info('remove by columns...')
    ex_filters = params['remove_by_colum']
    logging.info(pformat(ex_filters))
    for column, values in ex_filters.items():
        keep_indices.extend(
            adata[~adata.obs[column].isin(values)].obs_names
        )

# implicit filters
# mito filter
if 'mito_pct' in params:
    mito_threshold = params['mito_pct']
    logging.info(f'remove by mitochondrial percentage of {mito_threshold}%...')
    adata.var["mito"] = adata.var['feature_name'].str.startswith("MT-")
    adata.X = to_memory(adata.X)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    keep_indices.extend(
        adata[adata.obs['pct_counts_mito'] < mito_threshold].obs_names
    )

# remove sample with different number of cells per sample
if 'cells_per_sample' in params.keys():
    min_cells = params['cells_per_sample']['min']
    max_cells = params['cells_per_sample']['max']
    min_cells = min_cells if isinstance(min_cells, int) else 0
    max_cells = max_cells if isinstance(max_cells, int) else np.inf
    logging.info(f'remove by min={min_cells}, max={max_cells} cells per sample...')
    n_cells_per_sample = adata.obs['sample'].value_counts()
    samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)]
    keep_indices.extend(
        adata[adata.obs['sample'].isin(samples_to_keep)].obs_names
    )

keep_indices = set(keep_indices)
adata_removed = adata.n_obs - len(keep_indices)
logging.info(f'removed {adata_removed} cells')

if adata_removed == 0:
    save_empty(input_zarr, output_zarr, output_removed)

logging.info('Save filtered...')
adata[adata.obs_names.isin(keep_indices)].write_zarr(output_zarr)

logging.info('Save removed...')
adata[~adata.obs_names.isin(keep_indices)].write_zarr(output_removed)
