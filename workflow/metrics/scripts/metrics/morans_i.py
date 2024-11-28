import anndata as ad
import scanpy as sc
import numpy as np

from utils.accessors import subset_hvg

def get_morans_i_for_category_mean(_ad, _cov, num_mappings=10) -> float:
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
        morans.append(sc.metrics.morans_i(_ad, vals=mapped_values))
    
    # Calculate and return the mean of all Moran's I values
    return np.mean(morans)


# Moran's I for a random normal distribution (used as a negative baseline)
def morans_i_random(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    # Subset to highly variable genes
    adata.var[var_key] = adata.var[var_key]
    adata, _ = subset_hvg(adata, var_column=var_key, min_cells=1, compute_dask=True)

    try:
        return sc.metrics.morans_i(adata, vals=np.random.normal(size=adata.n_obs))
    except ZeroDivisionError:
        return float('nan')

# Moran's I for the batch (using permutations)
def morans_i_batch(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    # Subset to highly variable genes
    adata.var[var_key] = adata.var[var_key]
    adata, _ = subset_hvg(adata, var_column=var_key, min_cells=1, compute_dask=True)

    try:
        return get_morans_i_for_category_mean(adata, batch_key)
    except ZeroDivisionError:
        return float('nan')

#  Moran's I for the label/cell type (using permutations)
def morans_i_label(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    # Subset to highly variable genes
    adata.var[var_key] = adata.var[var_key]
    adata, _ = subset_hvg(adata, var_column=var_key, min_cells=1, compute_dask=True)

    try:
        return get_morans_i_for_category_mean(adata, label_key)
    except ZeroDivisionError:
        return float('nan')