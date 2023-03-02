import numpy as np
from .utils import select_neighbors


def isolated_label_f1(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    adata = select_neighbors(adata, output_type)

    return scib.me.isolated_labels_f1(
        adata,
        label_key=meta['label'],
        batch_key=meta['batch'],
        embed=None,
    )


def isolated_label_asw(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    adata = select_neighbors(adata, output_type)

    return scib.me.isolated_labels_asw(
        adata,
        label_key=meta['label'],
        batch_key=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )