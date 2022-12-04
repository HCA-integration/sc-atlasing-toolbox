import numpy as np
import pandas as pd
import scanpy as sc
import scib_metrics
from utils import write_metrics, get_from_adata, rename_categories

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

labels = rename_categories(adata, meta['label'])
batches = rename_categories(adata, meta['batch'])

output_types = meta['output_types']
scores = []
for output_type in output_types:
    score = np.nan
    if output_type == 'knn':
        continue

    if output_type == 'embed':
        X = adata.obsm['X_emb']
    else:
        X = adata.obsm['X_pca']

    X = X if isinstance(X, np.ndarray) else X.todense()

    score = scib_metrics.silhouette_batch(
        X=X,
        batch=batches,
        labels=labels,
    )
    scores.append(score)

write_metrics(scores, output_types, metric, metric_type, method, dataset, output_file)
