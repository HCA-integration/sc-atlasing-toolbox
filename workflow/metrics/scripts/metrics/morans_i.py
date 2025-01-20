import anndata as ad
import scanpy as sc
import numpy as np
from scipy import sparse
import pandas as pd
import itertools

from .bootstrap import bootstrap_metric
# from .utils_pipeline.processing import sc as rsc


def _morans_i(adata, covariate, **kwargs) -> float:
    """
    Moran's I wrapper
    Ensure the value is positive
    """
    if isinstance(covariate, str):
        values = adata.obs[covariate].values
    elif isinstance(covariate, pd.Series):
        values = covariate.values
    elif isinstance(covariate, (np.ndarray, sparse.spmatrix)):
        values = covariate.T
    else:
        raise ValueError(f'Invalid covariate type: {type(covariate)}')

    try:
        score = sc.metrics.morans_i(
            adata.obsp['connectivities'],
            vals=values,
            **kwargs
        )
    except ZeroDivisionError:
        return float('nan')
    
    if isinstance(score, np.ndarray):
        return [max(s, 0) for s in score]
    return max(score, 0)


def morans_i_categorical(adata, covariate, num_mappings=10) -> float:
    # Extract unique categories from the covariate column
    unique_categories = adata.obs[covariate].unique()
    num_categories = len(unique_categories)

    total_permutations = np.math.factorial(num_categories) # Determine if the number of mappings requested is feasible
    
    if num_mappings >= total_permutations:
        # If the requested number of mappings exceeds or equals the total possible permutations, generate all permutations
        unique_permutations = np.array(list(itertools.permutations(range(1, num_categories + 1))))
    else:
        # Otherwise, generate the required number of unique random permutations
        unique_permutations = np.array([np.random.permutation(range(1, num_categories + 1)) for _ in range(num_mappings)])

    morans = []  # Initialize a list to store Moran's I values
    for perm in unique_permutations:
        mapping_dict = dict(zip(unique_categories, perm))
        score = _morans_i(adata, covariate=adata.obs[covariate].map(mapping_dict))
        morans.append(score)
    
    # Calculate and return the mean of all Moran's I values
    return np.mean(morans)


def morans_i_random(adata, output_type, n_bootstraps=5, bootstrap_size=None, **kwargs):
    """
    Moran's I for a random normal distribution (used as a negative baseline)
    """
    return bootstrap_metric(
        adata,
        metric_function=lambda _ad: _morans_i(_ad, covariate=np.random.normal(size=_ad.n_obs)),
        n_bootstraps=n_bootstraps,
        size=bootstrap_size,
    )


def morans_i(adata, output_type, covariate, n_bootstraps=5, bootstrap_size=None, **kwargs):
    """
    Moran's I for covariate
    """
    covariates = covariate if isinstance(covariate, list) else [covariate]
    metrics_names = [f"M's I:{covariate}" for covariate in covariates]
    
    scores = []
    for covariate in covariates:
        is_numeric = pd.api.types.is_numeric_dtype(adata.obs[covariate])
        morans_i_func = _morans_i if is_numeric else morans_i_categorical
        score = morans_i_func(adata, covariate)
        # bootstrap_metric(
        #     adata,
        #     metric_function=morans_i_func,
        #     n_bootstraps=n_bootstraps,
        #     size=bootstrap_size,
        #     _cov=covariate,
        # )
        scores.append(score)
    
    return scores, metrics_names
