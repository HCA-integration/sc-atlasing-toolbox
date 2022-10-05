import numpy as np
import scanpy as sc
import scib
from utils import cluster_optimal, prepare_unintegrated

input_adata = snakemake.input.h5ad
output_scib = snakemake.output.metric

metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
label = snakemake.params.label
batch = snakemake.params.batch
metric_type = snakemake.params.metric_type
n_cores = snakemake.threads

print(f'read {input_adata} ...')
adata = sc.read(input_adata)
#adata.uns['log1p']["base"] = None  # quickfix, known bug

# evaluate only on labeled cells
adata = adata[adata.obs[label].notnull()]

# TODO: prepare integrated data depending on preprocessing metadata
# TODO: handle different output types

print(f'calculate {metric}...')
if metric == 'ari':
    res_max, score_max, score_all = cluster_optimal(
        adata,
        label,
        cluster_key='cluster',
        cluster_function=sc.tl.leiden,
        metric=scib.me.nmi,
        use_rep='X_emb',
        n_iterations=5
    )
    print(score_all)
    print(res_max, score_max)
    score = scib.me.ari(adata, 'cluster', label)
elif metric == 'asw_label':
    score = scib.me.silhouette(adata, label, embed='X_emb')
elif metric == 'asw_batch':
    score = scib.me.silhouette_batch(adata, batch, label, embed='X_emb', verbose=False)
elif metric == 'cell_cycle':
    adata_unintegrated = prepare_unintegrated(adata)
    adata_unintegrated.X = adata_unintegrated.X.todense()
    adata.X = adata.X.todense()
    score = scib.me.cell_cycle(adata_unintegrated, adata, batch, embed='X_emb', organism='human', verbose=False)
elif metric == 'clisi':
    score = scib.me.clisi_graph(
        adata,
        batch,
        label,
        type_='embed',
        subsample=50,
        n_cores=n_cores,
        verbose=False
    )
elif metric == 'graph_connectivity':
    score = scib.me.graph_connectivity(adata, label)
elif metric == 'ilisi':
    score = scib.me.ilisi_graph(
        adata,
        batch,
        type_='embed',
        subsample=50,
        n_cores=n_cores,
        verbose=False
    )
elif metric == 'isolated_label_asw':
    score = scib.me.isolated_labels(adata, label, batch, embed='X_emb', cluster=False, verbose=False)
elif metric == 'isolated_label_f1':
    score = scib.me.isolated_labels(adata, label, batch, embed='X_emb', cluster=True, verbose=False)
elif metric == 'nmi':
    #scib.me.opt_louvain(adata, label_key=label_key, cluster_key='cluster', use_rep='X_emb', verbose=False)
    res_max, score_max, score_all = cluster_optimal(
        adata,
        label,
        cluster_key='cluster',
        cluster_function=sc.tl.leiden,
        metric=scib.me.nmi,
        use_rep='X_emb',
        n_iterations=5
    )
    print(score_all)
    print(res_max, score_max)
    score = scib.me.nmi(adata, 'cluster', label)
elif metric == 'pcr':
    adata_unintegrated = prepare_unintegrated(adata)
    adata_unintegrated.X = adata_unintegrated.X.todense()
    adata.X = adata.X.todense()
    score = scib.me.pcr_comparison(adata_unintegrated, adata, batch, embed='X_emb', verbose=False)
else:
    print(f'Warning, {metric} is not defined. Returning NaN.')
    score = np.nan

print(f'{metric}: {score}')

# Save
with open(output_scib, 'w') as f:
    f.write(str(score))
