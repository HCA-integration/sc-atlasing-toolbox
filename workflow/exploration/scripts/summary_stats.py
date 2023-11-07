import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


input_zarr = snakemake.input.zarr
output = snakemake.output.tsv
output_sample = snakemake.output.sample
output_donor = snakemake.output.donor
study = snakemake.wildcards.study


logging.info(f'Read {input_zarr}...')
z = zarr.open(input_zarr)
obs = read_elem(z['obs'])
var = read_elem(z['var'])
n_obs = obs.shape[0]
n_vars = var.shape[0]

obs_per_sample = obs['sample'].value_counts()
if n_obs > 0:
    obs_per_sample.plot.barh(figsize=(10,8))
plt.savefig(output_sample)

obs_per_donor = obs['donor'].value_counts()
if n_obs > 0:
    obs_per_donor.plot.barh(figsize=(10,8))
plt.savefig(output_donor)

adata.var["mito"] = adata.var['feature_name'].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

stats = {
    'study': study,
    'n_cells': [n_obs],
    'n_genes': [n_vars],
    'n_samples': [obs_per_sample.shape[0]],
    'n_donors': [obs_per_donor.shape[0]],
    'median_cells_per_sample': [obs_per_sample.median()],
    'median_cells_per_donor': [obs_per_donor.median()],
    'disease_states': [','.join(obs['disease'].unique().to_list())],
    'mean_cells_per_sample': [obs['sample'].value_counts().mean()], 
    'median_counts_per_cell': [obs['total_counts'].median()],
    'median_genes_per_cell': [obs['n_genes_by_counts'].median()]
}

pd.DataFrame(stats).to_csv(output, sep='\t', index=False)
