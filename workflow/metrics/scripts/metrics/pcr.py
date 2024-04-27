import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from utils.accessors import adata_to_memory
from utils.assertions import assert_pca
from utils.misc import dask_compute

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


def pcr(adata, output_type, batch_key, label_key, adata_raw, n_threads=1, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    assert_pca(adata, check_varm=False)
    assert_pca(adata_raw, check_varm=False)
    # embed = 'X_emb' if output_type == 'embed' else 'X_pca'
    # assert embed in adata.obsm, f'Embedding {embed} missing from adata.obsm'
    
    return scib.me.pcr_comparison(
        adata_pre=adata_raw,
        adata_post=adata,
        covariate=batch_key,
        # embed=embed,  # assume that existing PCA is already computed on correct embedding
        n_threads=n_threads,
        linreg_method='numpy',
        recompute_pca=False,
    )


def pcr_y(adata, output_type, batch_key, label_key, adata_raw, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    adata_raw = adata_to_memory(adata_raw)
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

    # get correct feature names
    if 'feature_name' in adata_raw.var.columns:
        adata_raw.var_names = adata_raw.var['feature_name']

    # get organism
    upper_case_genes = sum(adata_raw.var_names.str.isupper())
    organism = 'mouse' if upper_case_genes <= 0.1 * adata_raw.n_vars else 'human'
    
    embed = 'X_emb' if output_type == 'embed' else 'X_pca'
    assert embed in adata.obsm, f'Embedding {embed} missing from adata.obsm'
    adata.obsm['X_emb'] = adata.obsm[embed]

    dask_compute(adata_raw)

    try:
        # compute cell cycle score per batch
        batch_key = batch_key
        for batch in adata.obs[batch_key].unique():
            scib.pp.score_cell_cycle(
                adata_raw[adata_raw.obs[batch_key] == batch],
                organism=organism
            )
    except Exception as e:
        raise ValueError(f'Error in score_cell_cycle: {e}') from e
    
    try:
        # compute score
        score = scib.me.cell_cycle(
            adata_pre=adata_raw,
            adata_post=adata,
            batch_key=batch_key,
            embed='X_emb',
            recompute_cc=False,
            organism=organism,
            verbose=False,
        )
    except ValueError as e:
        print(f'Warning: caught error in cell cycle score: {e}')
        score = np.nan

    return score
