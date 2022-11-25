"""
Data filtering
"""
import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params['filter']

adata = sc.read(input_h5ad)

# explicit filter
ex_filters = params['remove_by_colum']
for column, values in ex_filters.items():
    adata = adata[~adata.obs[column].isin(values)]

# implicit filters
# mito filter
if 'mito_pct' in params.keys():
    mito_threshold = params['mito_pct']
    adata.var["mito"] = adata.var['feature_name'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    adata = adata[adata.obs['pct_counts_mito'] < mito_threshold]

# remove sample with different number of cells per sample
if 'cells_per_sample' in params.keys():
    min_cells = params['cells_per_sample']['min']
    max_cells = params['cells_per_sample']['max']
    n_cells_per_sample = adata.obs['sample'].value_counts()
    samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)].index
    adata = adata[adata.obs['sample'].isin(samples_to_keep)]

adata.write(output_h5ad)
