"""
Data filtering
"""
from pprint import pformat
import numpy as np
import scanpy as sc
import anndata
import zarr
import logging
logging.basicConfig(level=logging.INFO)

from utils_pipeline.io import link_zarr


def save_empty(input_zarr, output_zarr, output_removed):
    # link unfiltered
    link_zarr(
        in_dir=input_zarr,
        out_dir=output_zarr,
        relative_path=False,
    )

    # save removed cells
    z = zarr.open(input_zarr)
    obs = anndata.experimental.read_elem(z['obs'])
    var = anndata.experimental.read_elem(z['var'])
    uns = anndata.experimental.read_elem(z['uns'])
    adata_empty = anndata.AnnData(
        X=np.zeros(shape=(0, var.shape[0])),
        obs=obs[:0],
        var=var,
        uns=uns
    )
    adata_empty.write_zarr(output_removed)
    # link_zarr(
    #     in_dir=input_zarr,
    #     out_dir=output_removed,
    #     file_names=['uns'],
    #     overwrite=True,
    #     relative_path=False,
    # )


input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
output_removed = snakemake.output.removed
params = snakemake.params.get('filter', {})

logging.info(pformat(params))

if len(params) == 0:
    save_empty(input_zarr, output_zarr, output_removed)
    exit(0)

adata = anndata.read_zarr(input_zarr)
logging.info(adata.__str__())

# subset to available annotations TODO: make this more generic
if 'cell_type' in adata.obs.columns:
    adata_filtered = adata[adata.obs['cell_type'].notna()]

# explicit filter
if 'remove_by_colum' in params:
    logging.info('remove by columns...')
    ex_filters = params['remove_by_colum']
    logging.info(pformat(ex_filters))
    for column, values in ex_filters.items():
        adata_filtered = adata_filtered[~adata_filtered.obs[column].isin(values)]

# implicit filters
# mito filter
if 'mito_pct' in params.keys():
    mito_threshold = params['mito_pct']
    logging.info(f'remove by mitochondrial percentage of {mito_threshold}%...')
    adata_filtered.var["mito"] = adata_filtered.var['feature_name'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata_filtered, qc_vars=["mito"], inplace=True)
    adata_filtered = adata_filtered[adata_filtered.obs['pct_counts_mito'] < mito_threshold]

# remove sample with different number of cells per sample
if 'cells_per_sample' in params.keys():
    min_cells = params['cells_per_sample']['min']
    max_cells = params['cells_per_sample']['max']
    min_cells = min_cells if isinstance(min_cells, int) else 0
    max_cells = max_cells if isinstance(max_cells, int) else np.inf
    logging.info(f'remove by min={min_cells}, max={max_cells} cells per sample...')
    n_cells_per_sample = adata_filtered.obs['sample'].value_counts()
    samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)].index
    adata_filtered = adata_filtered[adata_filtered.obs['sample'].isin(samples_to_keep)]

adata_removed = adata.n_obs - sum(adata.obs_names.isin(adata_filtered.obs_names))
logging.info(f'removed {adata_removed} cells')

if adata_removed == 0:
    save_empty(input_zarr, output_zarr, output_removed)
    exit(0)

logging.info('Save filtered...')
logging.info(adata_filtered.__str__())
adata_filtered.write_zarr(output_zarr)

logging.info('Save removed...')
adata[~adata.obs_names.isin(adata_filtered.obs_names)].write_zarr(output_removed)
