"""
Data filtering
"""
from pprint import pprint
import numpy as np
import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
output_removed = snakemake.output.removed
params = snakemake.params['filter']

adata = sc.read(input_h5ad)

# subset to available annotations
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
    min_cells = 0 if min_cells is None else min_cells
    max_cells = np.inf if max_cells is None else max_cells
    print(f'remove by min={min_cells}, max={max_cells} cells per sample...')
    n_cells_per_sample = adata_filtered.obs['sample'].value_counts()
    samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)].index
    adata_filtered = adata_filtered[adata_filtered.obs['sample'].isin(samples_to_keep)]

print('Save filtered...')
adata_filtered.write(output_h5ad)

adata_removed = adata[~adata.obs_names.isin(adata_filtered.obs_names)]
print(f'removed {adata_removed.n_obs} cells')

print('Save removed...')
adata.write(output_removed)
