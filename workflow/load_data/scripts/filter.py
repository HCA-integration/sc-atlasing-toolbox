"""
Data filtering
"""
from pprint import pprint
import numpy as np
import scanpy as sc
import anndata

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
output_removed = snakemake.output.removed
params = snakemake.params['filter']

adata = anndata.read_zarr(input_zarr)
print(adata)

# subset to available annotations TODO: make this more generic
if 'cell_type' in adata.obs.columns:
    adata_filtered = adata[adata.obs['cell_type'].notna()]

# explicit filter
if 'remove_by_colum' in params:
    print('remove by columns...')
    ex_filters = params['remove_by_colum']
    pprint(ex_filters)
    for column, values in ex_filters.items():
        adata_filtered = adata_filtered[~adata_filtered.obs[column].isin(values)]

# implicit filters
# mito filter
if 'mito_pct' in params.keys():
    mito_threshold = params['mito_pct']
    print(f'remove by mitochondrial percentage of {mito_threshold}%...')
    adata_filtered.var["mito"] = adata_filtered.var['feature_name'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata_filtered, qc_vars=["mito"], inplace=True)
    adata_filtered = adata_filtered[adata_filtered.obs['pct_counts_mito'] < mito_threshold]

# remove sample with different number of cells per sample
if 'cells_per_sample' in params.keys():
    min_cells = params['cells_per_sample']['min']
    max_cells = params['cells_per_sample']['max']
    min_cells = min_cells if isinstance(min_cells, int) else 0
    max_cells = max_cells if isinstance(max_cells, int) else np.inf
    print(f'remove by min={min_cells}, max={max_cells} cells per sample...')
    n_cells_per_sample = adata_filtered.obs['sample'].value_counts()
    samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)].index
    adata_filtered = adata_filtered[adata_filtered.obs['sample'].isin(samples_to_keep)]

print('Save filtered...')
print(adata_filtered)
adata_filtered.write_zarr(output_zarr)

adata_removed = adata.n_obs - sum(adata.obs_names.isin(adata_filtered.obs_names))
print(f'removed {adata_removed} cells')

print('Save removed...')
adata[~adata.obs_names.isin(adata_filtered.obs_names)].write_zarr(output_removed)
