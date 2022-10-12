import pandas as pd
import scanpy as sc
import scib
from utils import write_metrics, get_from_adata, compute_neighbors

input_adata = snakemake.input.h5ad
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
n_cores = snakemake.threads

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
meta = get_from_adata(adata)

records = []
for output_type in meta['output_types']:
    adata = compute_neighbors(adata, output_type)
    score = scib.me.isolated_labels(
        adata,
        label_key=meta['label'],
        batch_key=meta['batch'],
        embed=None,
        cluster=True,
    )
    records.append((metric, method, output_type, metric_type, score))

write_metrics(records, output_file)
