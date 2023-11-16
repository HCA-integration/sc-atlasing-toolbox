import pandas as pd
import scanpy as sc

from utils.io import read_anndata

input_zarr = snakemake.input.zarr
output_obs = snakemake.output.obs

adata = read_anndata(snakemake.input[0])

print('Calculate QC stats...')
if 'feature_name' in adata.var.columns:
    var_names = adata.var['feature_name']
else:
    var_names = adata.var_names

adata.var["mito"] = var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

adata.obs.to_csv(output_obs, sep='\t')
