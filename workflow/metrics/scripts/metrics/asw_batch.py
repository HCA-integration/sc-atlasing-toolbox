import numpy as np
import pandas as pd
import scanpy as sc
import scib
from utils import write_metrics, get_from_adata

input_adata = snakemake.input.h5ad
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
dataset = snakemake.wildcards.dataset
hyperparams = snakemake.wildcards.hyperparams

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
meta = get_from_adata(adata)

# evaluate only on labeled cells
# adata = adata[adata.obs[label].notnull()]

output_types = meta['output_types']
scores = []
for output_type in output_types:
    score = np.nan
    if output_type == 'knn':
        continue

    score = scib.me.silhouette_batch(
        adata,
        batch_key=meta['batch'],
        label_key=meta['label'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )
    scores.append(score)

write_metrics(
    scores=scores,
    output_types=output_types,
    metric=metric,
    metric_type=metric_type,
    method=method,
    hyperparams=hyperparams,
    dataset=dataset,
    filename=output_file
)
