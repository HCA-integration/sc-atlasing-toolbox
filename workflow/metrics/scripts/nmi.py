import pandas as pd
import scanpy as sc
import scib
from utils import cluster_optimal

input_adata = snakemake.input.h5ad
output_scib = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
label = adata.uns['integration']['label_key']
output_type = adata.uns['integration']['output_type']
output_types = [output_type] if isinstance(output_type, str) else output_type

# evaluate only on labeled cells
# adata = adata[adata.obs[label].notnull()]

rep_map = {
    'knn': None,
    'embed': 'X_emb',
    'full': 'X_pca'
}

records = []
for output_type in output_types:
    _, score, _ = cluster_optimal(
        adata,
        label,
        cluster_key='cluster',
        cluster_function=sc.tl.leiden,
        metric=scib.me.nmi,
        use_rep=rep_map[output_type],
        n_iterations=5
    )
    records.append((metric, method, output_type, metric_type, score))

df = pd.DataFrame.from_records(
    records,
    columns=['metric', 'method', 'output_type', 'metric_type', 'score']
)
df.to_csv(output_scib, sep='\t', index=False)
