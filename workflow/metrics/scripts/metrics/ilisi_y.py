import pandas as pd
import scanpy as sc
import scib_metrics
from utils import write_metrics, get_from_adata, rename_categories, compute_neighbors

input_adata = snakemake.input.h5ad
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
dataset = snakemake.wildcards.dataset
hyperparams = snakemake.wildcards.hyperparams
n_cores = snakemake.threads

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
meta = get_from_adata(adata)
adata_raw = adata.raw.to_adata()

batches = rename_categories(adata, meta['batch'])

output_types = meta['output_types']
scores = []
for output_type in output_types:
    adata = compute_neighbors(adata, output_type)
    score = scib_metrics.ilisi_knn(
        X=adata.obsp['distances'].toarray(),
        batches=batches
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
