import scanpy as sc
from .utils import cluster_optimal, select_neighbors, rename_categories


def ari(adata, output_type, meta):
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
    return scib.me.ari(adata, meta['label'], 'cluster')


def ari_leiden_y(adata, output_type, meta):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, meta['label'])

    _, score = scib_metrics.nmi_ari_cluster_labels_leiden(
        X=adata.obsp['connectivities'].toarray(),
        labels=labels,
        optimize_resolution=True,
    )
    return score


def ari_kmeans_y(adata, output_type, meta):
    import scib_metrics

    labels = rename_categories(adata, meta['label'])
    adata = select_neighbors(adata, output_type)
    _, score = scib_metrics.nmi_ari_cluster_labels_kmeans(
        X=adata.obsp['connectivities'].toarray(),
        labels=labels,
    )
    return score
