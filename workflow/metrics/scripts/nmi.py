import pandas as pd
import scanpy as sc
import scib
from utils import cluster_optimal

input_adata = snakemake.input.h5ad
output_scib = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
params = snakemake.params

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
label = adata.uns['integration']['label_key']
ot = adata.uns['integration']['output_type']
output_types = [ot] if isinstance(ot, str) else ot

# evaluate only on labeled cells
#adata = adata[adata.obs[label].notnull()]

for ot in output_types:
    if ot in ['full', 'embed']:
        assert 'X_emb' in adata.obsm
        # cluster on graph built from 'X_emb'
        res_max, score_max, score_all = cluster_optimal(
            adata,
            label,
            cluster_key='cluster',
            cluster_function=sc.tl.leiden,
            metric=scib.me.nmi,
            use_rep='X_emb',
            n_iterations=5
        )
    elif ot == 'knn':
        # clustering on the knn graph directly
        res_max, score_max, score_all = cluster_optimal(
            adata,
            label,
            cluster_key='cluster',
            cluster_function=sc.tl.leiden,
            metric=scib.me.nmi,
            use_rep=None,
            n_iterations=5
        )
    else:
        raise ValueError(f'invalid output type: {ot}')

    score = scib.me.nmi(adata, 'cluster', label)

    df = pd.DataFrame(
        {
            'metric': metric,
            'method': method,
            'output_type': ot,
            'metric_type': params['metric_type'],
            'score': score
        }
    )
    df.to_csv(output_scib, sep='\t')
