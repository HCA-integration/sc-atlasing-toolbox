"""
Highly variable gene selection
- lineage specific HVGs
"""
import anndata
import scanpy as sc

input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
n_hvgs = snakemake.params['n_hvgs']
batch_key = snakemake.params['batch']

adata = anndata.read_zarr(input_file)
adata.uns["log1p"] = {"base": None}

# TODO: lineage specific (lineage + batch)
adata.obs['hvg_batch'] = adata.obs[batch_key]
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, batch_key='hvg_batch')

adata.write_zarr(output_file)
