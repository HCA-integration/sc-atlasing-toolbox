import numpy as np
import anndata as ad
import scanpy as sc
from .utils_pipeline.processing import sc as rsc


def get_bootstrap_adata(adata: ad.AnnData, size: int = None, use_rep='X_emb'):
    if size is None:
        size = int(0.9 * adata.n_obs)
    
    print(f'Subset to {size} cells...', flush=True)
    rand_indices = np.random.choice(adata.obs_names, size=size, replace=False)
    adata_bootstrap = adata[rand_indices].copy()
    
    # Recalculation is important as kNN graph is input for Moran's I
    print('Recalculate neighbors from distances...', flush=True)
    n_neighbors = adata.uns.get('neighbors', {}).get('params', {}).get('n_neighbors', 15)
    try:
        adata_bootstrap.obsm['distances'] = adata_bootstrap.obsp['distances']
        sc.pp.neighbors(
            adata_bootstrap,
            use_rep='distances',
            metric='precomputed',
            transformer='sklearn',
            n_neighbors=n_neighbors,
        )
    except ValueError as e:
        print(f'Error in recalculation of neighbors:\n{e}', flush=True)
        print('Recalculate neighbors from {use_rep}...', flush=True)
        rsc.pp.neighbors(
            adata_bootstrap,
            use_rep=use_rep,
            n_neighbors=n_neighbors
        )
    finally:
        return adata_bootstrap


def bootstrap_metric(
    adata: ad.AnnData,
    metric_function,
    n_bootstraps: int = 5,
    size: int = None,
    subset_kwargs: tuple = (),
    *args,
    **kwargs,
):
    """
    Call a metric function for multiple bootstraps
    
    :param n_bootstraps: number of bootstraps if multiple subsetting is desired
    :return: single score if n_bootstraps <= 1, else list of scores of length n_bootstraps
    """
    
    if n_bootstraps <= 1:
        return metric_function(adata, *args, **kwargs)
    
    scores = []
    for _ in range(n_bootstraps):
        _ad = get_bootstrap_adata(adata, size=size)
        score = metric_function(_ad, *args, **kwargs)
        scores.append(score)
    
    return scores