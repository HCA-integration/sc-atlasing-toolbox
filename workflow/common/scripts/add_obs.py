"""
Add columns to adata.obs
"""
import pandas as pd
import scanpy as sc

input_h5ad = snakemake.input.h5ad
input_tsv = snakemake.input.tsv
output_h5ad = snakemake.output.h5ad

tsv_column = snakemake.params['tsv_column']
h5ad_column = snakemake.params['h5ad_column']


annotation = pd.read_table(input_tsv)
adata = sc.read(input_h5ad)

annotation.index = annotation[tsv_column]
adata.obs=pd.merge(adata.obs, annotation)

adata.write(output_h5ad)
