import numpy as np
from .utils import rename_categories


def asw_batch(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    return scib.me.silhouette_batch(
        adata,
        batch_key=meta['batch'],
        group_key=meta['label'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def asw_batch_y(adata, output_type, meta):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    if output_type == 'embed':
        X = adata.obsm['X_emb']
    else:
        X = adata.obsm['X_pca']

    X = X if isinstance(X, np.ndarray) else X.todense()
    labels = rename_categories(adata, meta['label'])
    batches = rename_categories(adata, meta['batch'])

    return scib_metrics.silhouette_batch(
        X=X,
        batch=batches,
        labels=labels,
    )


def asw_label(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    return scib.me.silhouette(
        adata,
        group_key=meta['label'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def asw_label_y(adata, output_type, meta):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    if output_type == 'embed':
        X = adata.obsm['X_emb']
    else:
        X = adata.obsm['X_pca']

    X = X if isinstance(X, np.ndarray) else X.todense()
    labels = rename_categories(adata, meta['label'])

    return scib_metrics.silhouette_label(
        X=X,
        labels=labels,
    )
