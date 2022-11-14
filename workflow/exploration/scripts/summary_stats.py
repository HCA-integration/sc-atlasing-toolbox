import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt

input_h5ad = snakemake.input.h5ad
output = snakemake.output.tsv
output_sample = snakemake.output.sample
output_donor = snakemake.output.donor

print(f'Read {input_h5ad}...')
adata = sc.read(input_h5ad)

obs_per_sample = adata.obs['sample'].value_counts()
obs_per_sample.plot.barh(figsize=(10,8))
plt.savefig(output_sample)

obs_per_donor = adata.obs['donor'].value_counts()
obs_per_donor.plot.barh(figsize=(10,8))
plt.savefig(output_donor)

stats = {
    'dataset': snakemake.wildcards.dataset,
    'n_cells': [adata.n_obs],
    'n_genes': [adata.n_vars],
    'n_samples': [obs_per_sample.shape[0]],
    'n_donors': [obs_per_donor.shape[0]],
    'median_cells_per_sample': [obs_per_sample.median()],
    'median_cells_per_donor': [obs_per_donor.median()],
    'disease_states': [','.join(adata.obs['disease'].unique().to_list())],
}

pd.DataFrame(stats).to_csv(output, sep='\t', index=False)
