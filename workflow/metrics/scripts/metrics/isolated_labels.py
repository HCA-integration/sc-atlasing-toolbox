import numpy as np
from .utils import select_neighbors


def isolated_label_f1(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    return scib.me.isolated_labels_f1(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed=None,
    )


def isolated_label_asw(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    return scib.me.isolated_labels_asw(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def isolated_label_asw_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else X.todense()

    return scib_metrics.isolated_labels(
        X=X,
        labels=adata.obs[label_key],
        batch=adata.obs[batch_key],
    )
