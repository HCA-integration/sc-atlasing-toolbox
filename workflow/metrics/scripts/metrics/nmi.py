import numpy as np
import scanpy as sc
from .utils import cluster_optimal, select_neighbors, rename_categories


def nmi(adata, output_type, meta):
    import scib

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
    return scib.me.nmi(adata, meta['label'], 'cluster')


def nmi_leiden_y(adata, output_type, meta):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, meta['label'])

    scores = scib_metrics.nmi_ari_cluster_labels_leiden(
        X=adata.obsp['connectivities'],
        labels=labels,
        optimize_resolution=True,
    )
    return scores['nmi']


def nmi_kmeans_y(adata, output_type, meta):
    import scib_metrics

    labels = rename_categories(adata, meta['label'])
    adata = select_neighbors(adata, output_type)

    X = adata.obsp['connectivities']
    X = X if isinstance(X, np.ndarray) else X.to_array()

    scores = scib_metrics.nmi_ari_cluster_labels_kmeans(
        X=X,
        labels=labels,
    )
    return scores['nmi']
