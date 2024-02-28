import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
output_sample = snakemake.output.sample
output_donor = snakemake.output.donor
study = snakemake.wildcards.file_id
sample = snakemake.params.sample
donor = snakemake.params.donor
categories = snakemake.params.get('categories', [])
if isinstance(categories, str):
    categories = [categories]


logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', var='var')

obs_per_sample = adata.obs[sample].value_counts()
if adata.n_obs > 0:
    obs_per_sample.plot.barh(figsize=(10,8))
plt.savefig(output_sample)

obs_per_donor = adata.obs[donor].value_counts()
if adata.n_obs > 0:
    obs_per_donor.plot.barh(figsize=(10,8))
plt.savefig(output_donor)

adata.var["mito"] = adata.var['feature_name'].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

stats = {
    'study': study,
    'n_cells': [adata.n_obs],
    'n_genes': [adata.n_vars],
    'n_samples': [obs_per_sample.shape[0]],
    'n_donors': [obs_per_donor.shape[0]],
    'median_cells_per_sample': [obs_per_sample.median()],
    'median_cells_per_donor': [obs_per_donor.median()],
    'disease_states': [','.join(obs['disease'].unique().to_list())],
    'mean_cells_per_sample': [obs['sample'].value_counts().mean()], 
    'median_counts_per_cell': [obs['total_counts'].median()],
    'median_genes_per_cell': [obs['n_genes_by_counts'].median()],
}

for cat in categories:
    assert cat in adata.obs.columns, f'{cat} not in {adata.obs.columns}'
    stats |= {cat: [','.join(adata.obs[cat].value_counts().index.tolist())]}

pd.DataFrame(stats).to_csv(output_tsv, sep='\t', index=False)
