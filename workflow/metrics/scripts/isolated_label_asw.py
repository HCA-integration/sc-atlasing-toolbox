import numpy as np
import pandas as pd
import scanpy as sc
import scib
from utils import write_metrics, get_from_adata, compute_neighbors

input_adata = snakemake.input.h5ad
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
dataset = snakemake.wildcards.dataset

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
meta = get_from_adata(adata)

output_types = meta['output_types']
scores = []
for output_type in output_types:
    score = np.nan
    if output_type == 'knn':
        continue

    score = scib.me.isolated_labels(
        adata,
        label_key=meta['label'],
        batch_key=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
        cluster=False,
    )
    scores.append(score)

write_metrics(scores, output_types, metric, metric_type, method, dataset, output_file)