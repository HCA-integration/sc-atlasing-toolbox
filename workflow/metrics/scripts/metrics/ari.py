import pandas as pd
import scanpy as sc
import scib
from utils import cluster_optimal, write_metrics, get_from_adata, select_neighbors

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

# # deprecated bs using precomputed neighbors
# rep_map = {
#     'knn': None,
#     'embed': 'X_emb',
#     'full': 'X_pca'
# }

output_types = meta['output_types']
scores = []
for output_type in output_types:
    adata = select_neighbors(adata, output_type)
    cluster_optimal(
        adata=adata,
        label_key=meta['label'],
        cluster_key='cluster',
        cluster_function=sc.tl.leiden,
        metric=scib.me.nmi,
        use_rep=None,
        n_iterations=5
    )
    score = scib.me.ari(adata, meta['label'], 'cluster')
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
