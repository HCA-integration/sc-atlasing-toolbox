import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from utils.assertions import assert_pca
from utils.misc import dask_compute

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


def pcr_comparison(adata, output_type, batch_key, label_key, adata_raw, n_threads=1, **kwargs):
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
        # n_threads=n_threads,
        linreg_method='numpy',
        recompute_pca=False,
    )


def pcr_y(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    adata_raw = dask_compute(adata_raw[:, adata_raw.var[var_key]].copy())
    X_pre = adata_raw.X
    X_post = adata.obsm['X_emb'] if output_type == 'embed' else adata.X
    X_pre, X_post = [X if isinstance(X, np.ndarray) else X.todense() for X in [X_pre, X_post]]

    return scib_metrics.pcr_comparison(
        X_pre=X_pre,
        X_post=X_post,
        covariate=adata.obs[batch_key],
        categorical=True
    )


def cell_cycle(
    adata,
    output_type,
    batch_key,
    label_key,
    adata_raw,
    var_key='metrics_features',
    n_threads=1,
    **kwargs
):
    import scib

    if output_type == 'knn':
        return np.nan

    # subset adata if dataset too large
    n_subset = int(1e6)
    if adata.n_obs > n_subset:
        random_subset = np.random.choice(adata.obs_names, size=n_subset, replace=False)
        adata = adata[random_subset].copy()
        adata_raw = adata_raw[random_subset].copy()

    # get correct feature names
    if 'feature_name' in adata_raw.var.columns:
        adata_raw.var_names = adata_raw.var['feature_name']

    # get organism
    upper_case_genes = sum(adata_raw.var_names.str.isupper())
    organism = 'mouse' if upper_case_genes <= 0.1 * adata_raw.n_vars else 'human'
    
    embed = 'X_emb' if output_type == 'embed' else 'X_pca'
    assert embed in adata.obsm, f'Embedding {embed} missing from adata.obsm'
    adata.obsm['X_emb'] = adata.obsm[embed]

    # TODO: ensure that cell cycle genes are preserved
    adata_raw = dask_compute(adata_raw[:, adata_raw.var[var_key]].copy())

    # compute cell cycle score per batch
    for batch in adata.obs[batch_key].unique():
        batch_mask = adata_raw.obs[batch_key] == batch
        try:
            ad_sub = adata_raw[batch_mask].copy()
            scib.pp.score_cell_cycle(ad_sub, organism=organism)
            adata_raw.obs.loc[batch_mask, 'S_score'] = ad_sub.obs['S_score']
            adata_raw.obs.loc[batch_mask, 'G2M_score'] = ad_sub.obs['G2M_score']
        except Exception as e:
            print(f'Error in score_cell_cycle for batch={batch}:\n{e}', flush=True)
            adata_raw.obs.loc[batch_mask, 'S_score'] = 0
            adata_raw.obs.loc[batch_mask, 'G2M_score'] = 0
    
    # try:
    # compute score
    score = scib.me.cell_cycle(
        adata_pre=adata_raw,
        adata_post=adata,
        batch_key=batch_key,
        embed='X_emb',
        recompute_cc=False,
        organism=organism,
        verbose=False,
        linreg_method='sklearn',
        n_threads=n_threads,
    )
    # except ValueError as e:
    #     print(f'Warning: caught error in cell cycle score: {e}')
    #     score = np.nan

    return score


def pcr_random(adata, output_type, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    assert_pca(adata, check_varm=False)

    # Add column full of random values and use it as input for the PCR
    adata.obs['random'] = np.random.normal(size=adata.n_obs)

    return scib.metrics.pcr(
        adata,
        covariate='random',
        recompute_pca=False,
        linreg_method='numpy',
    )
    

def pcr(adata, output_type, covariate, **kwargs):
    import scib
    
    covariates = covariate if isinstance(covariate, list) else [covariate]
    metric_names = [f'pcr:{cov}' for cov in covariates]

    if output_type == 'knn':
        return [np.nan] * len(covariates), metric_names, 

    assert_pca(adata, check_varm=False)

    scores = [
        scib.metrics.pcr(
            adata,
            covariate=covariate,
            recompute_pca=False,
            linreg_method='numpy',
        )
        for covariate in covariates
    ]
    
    return scores, metric_names
