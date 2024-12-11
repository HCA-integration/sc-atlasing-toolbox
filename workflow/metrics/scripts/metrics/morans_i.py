import anndata as ad
import scanpy as sc
import numpy as np

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


def get_bootstrap_adata(adata, size):
    rand_indices = np.random.choice(adata.obs.index, size=size, replace=False)
    adata_bootstrap = adata[adata.obs.index.isin(rand_indices)]
    sc.pp.neighbors(adata_bootstrap, use_rep='X_pca')   # Recalculation is important as kNN graph is input for Moran's I
    return adata_bootstrap


# Moran's I for a random normal distribution (used as a negative baseline)
def morans_i_random(adata, output_type, **kwargs):
    morans_values = []

    for _ in range(5):
        ad = get_bootstrap_adata(adata, int(0.9*adata.n_obs))
        try:
            value = sc.metrics.morans_i(ad, vals=np.random.normal(size=ad.n_obs))
        except ZeroDivisionError:
            value = float('nan')

        morans_values.append(value)
    
    return morans_values


# Moran's I for the batch (using permutations)
def morans_i_batch(adata, output_type, batch_key, **kwargs):
    morans_values = []

    for _ in range(5):
        ad = get_bootstrap_adata(adata, int(0.9*adata.n_obs))
        morans_values.append(get_morans_i_for_category_mean(ad, batch_key))
    
    return morans_values


#  Moran's I for the label/cell type (using permutations)
def morans_i_label(adata, output_type, label_key, **kwargs):
    morans_values = []

    for _ in range(5):
        ad = get_bootstrap_adata(adata, int(0.9*adata.n_obs))
        morans_values.append(get_morans_i_for_category_mean(ad, label_key))
    
    return morans_values