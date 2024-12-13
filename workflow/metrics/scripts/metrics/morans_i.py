import anndata as ad
import scanpy as sc
import numpy as np

from .bootstrap import bootstrap_metric


def morans_i(*args, **kwargs):
    """
    Moran's I wrapper
    Ensure the value is positive
    """
    # try:
    score = sc.metrics.morans_i(*args, **kwargs)
    return max(score, 0)
    # except ZeroDivisionError:
    #     return float('nan')


def morans_i_categorical(_ad, _cov, num_mappings=10) -> float:
    import itertools

    # Extract unique categories from the covariate column
    unique_categories = _ad.obs[_cov].unique()
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
        mapped_values = _ad.obs[_cov].map(mapping_dict).values
        score = morans_i(_ad, vals=mapped_values)
        morans.append(score)
    
    # Calculate and return the mean of all Moran's I values
    return np.mean(morans)


def morans_i_random(adata, output_type, n_bootstraps=5, bootstrap_size=None, **kwargs):
    """
    Moran's I for a random normal distribution (used as a negative baseline)
    """
    return bootstrap_metric(
        adata,
        metric_function=lambda _ad: morans_i(_ad, vals=np.random.normal(size=_ad.n_obs)),
        n_bootstraps=n_bootstraps,
        size=bootstrap_size,
    )


def morans_i_batch(adata, output_type, batch_key, n_bootstraps=5, bootstrap_size=None, **kwargs):
    """
    Moran's I for the batch (using permutations)
    """
    return bootstrap_metric(
        adata,
        metric_function=morans_i_categorical,
        n_bootstraps=n_bootstraps,
        size=bootstrap_size,
        _cov=batch_key,
    )


def morans_i_label(adata, output_type, label_key, n_bootstraps=5, bootstrap_size=None, **kwargs):
    """
    Moran's I for the label/cell type (using permutations)
    """
    return bootstrap_metric(
        adata,
        metric_function=morans_i_categorical,
        n_bootstraps=n_bootstraps,
        size=bootstrap_size,
        _cov=label_key,
    )
