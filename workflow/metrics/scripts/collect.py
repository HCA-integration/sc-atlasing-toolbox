import pandas as pd

from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
metrics_file = snakemake.input.metrics
output_file = snakemake.output[0]

adata = read_anndata(input_file, uns='uns')
metrics_df = pd.read_table(metrics_file, sep='\t')
for col in metrics_df.select_dtypes(include=['object', 'category']).columns:
    metrics_df[col] = metrics_df[col].astype(str)
adata.uns['metrics'] = metrics_df

write_zarr_linked(adata, input_file, output_file, files_to_keep=['uns'])
