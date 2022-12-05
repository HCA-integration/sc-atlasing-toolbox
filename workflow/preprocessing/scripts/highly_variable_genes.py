"""
Highly variable gene selection
- lineage specific HVGs
"""
import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
n_hvgs = snakemake.params['n_hvgs']
batch_key = snakemake.params['batch']

adata = sc.read(input_h5ad)
adata.uns["log1p"] = {"base": None}

# TODO: lineage specific (lineage + batch)
adata.obs['hvg_batch'] = adata.obs[batch_key]
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, batch_key='hvg_batch')

adata.write(output_h5ad)
