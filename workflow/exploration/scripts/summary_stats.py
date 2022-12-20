import anndata
import pandas as pd
from matplotlib import pyplot as plt

input_zarr = snakemake.input.zarr
output = snakemake.output.tsv
output_sample = snakemake.output.sample
output_donor = snakemake.output.donor
study = snakemake.wildcards.study


print(f'Read {input_zarr}...')
adata = anndata.read_zarr(input_zarr)

obs_per_sample = adata.obs['sample'].value_counts()
if adata.n_obs > 0:
    obs_per_sample.plot.barh(figsize=(10,8))
plt.savefig(output_sample)

obs_per_donor = adata.obs['donor'].value_counts()
if adata.n_obs > 0:
    obs_per_donor.plot.barh(figsize=(10,8))
plt.savefig(output_donor)

stats = {
    'study': study,
    'n_cells': [adata.n_obs],
    'n_genes': [adata.n_vars],
    'n_samples': [obs_per_sample.shape[0]],
    'n_donors': [obs_per_donor.shape[0]],
    'median_cells_per_sample': [obs_per_sample.median()],
    'median_cells_per_donor': [obs_per_donor.median()],
    'disease_states': [','.join(adata.obs['disease'].unique().to_list())],
}

pd.DataFrame(stats).to_csv(output, sep='\t', index=False)
