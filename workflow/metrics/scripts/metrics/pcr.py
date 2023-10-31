import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse


def pcr(adata, output_type, batch_key, label_key, adata_raw, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    return scib.me.pcr_comparison(
        adata_pre=adata_raw,
        adata_post=adata,
        covariate=batch_key,
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def pcr_y(adata, output_type, batch_key, label_key, adata_raw, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    X_pre = adata_raw.X

    X_post = adata.obsm['X_emb'] if output_type == 'embed' else adata.X
    X_pre, X_post = [X if isinstance(X, np.ndarray) else X.todense() for X in [X_pre, X_post]]

    return scib_metrics.pcr_comparison(
        X_pre=X_pre,
        X_post=X_post,
        covariate=adata.obs[batch_key],
        categorical=True
    )


def cell_cycle(adata, output_type, batch_key, label_key, adata_raw, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    if 'feature_name' in adata_raw.var.columns:
        adata_raw.var_names = adata_raw.var['feature_name']
    upper_case_genes = sum(adata.var_names.str.isupper())
    organism = 'mouse' if upper_case_genes <= 0.1 * adata.n_vars else 'human'

    try:
        # compute cell cycle score per batch
        batch_key = batch_key
        for batch in adata.obs[batch_key].unique():
            scib.pp.score_cell_cycle(
                adata_raw[adata_raw.obs[batch_key] == batch],
                organism=organism
            )
        
        # compute score
        score = scib.me.cell_cycle(
            adata_pre=adata_raw,
            adata_post=adata,
            batch_key=batch_key,
            embed='X_emb' if output_type == 'embed' else 'X_pca',
            recompute_cc=False,
            organism=organism,
            verbose=False,
        )
    except ValueError as e:
        warnings.warn(f'ValueError in cell cycle score: {e}')
        score = np.nan

    return score
