import pandas as pd
import scanpy as sc
import scib
from utils import cluster_optimal, write_metrics, get_from_adata

input_adata = snakemake.input.h5ad
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
meta = get_from_adata(adata)

# evaluate only on labeled cells
# adata = adata[adata.obs[label].notnull()]

rep_map = {
    'knn': None,
    'embed': 'X_emb',
    'full': 'X_pca'
}

records = []
for output_type in meta['output_types']:
    cluster_optimal(
        adata=adata,
        label_key=meta['label'],
        cluster_key='cluster',
        cluster_function=sc.tl.leiden,
        metric=scib.me.nmi,
        use_rep=rep_map[output_type],
        n_iterations=5
    )
    score = scib.me.ari(adata, meta['label'], 'cluster')
    records.append((metric, method, output_type, metric_type, score))

write_metrics(records, output_file)
