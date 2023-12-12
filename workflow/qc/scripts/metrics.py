import pandas as pd
import scanpy as sc
import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, link_zarr_partial

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

adata = read_anndata(snakemake.input[0], X='X', obs='obs', var='var')

if adata.n_obs == 0:
    logging.info(f'Write empty zarr file to {output_zarr}...')
    ad.AnnData(obs=adata.obs).write_zarr(output_zarr)
    exit(0)

print('Calculate QC stats...')
if 'feature_name' in adata.var.columns:
    var_names = adata.var['feature_name']
else:
    var_names = adata.var_names

adata.var["mito"] = var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

logging.info(f'Write zarr file to {output_zarr}...')
ad.AnnData(obs=adata.obs).write_zarr(output_zarr)
link_zarr_partial(input_zarr, output_zarr, files_to_keep=['obs'])