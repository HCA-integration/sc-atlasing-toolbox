import anndata as ad
import scib 
import numpy as np

def pcr_random(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Add column full of random values and use it as input for the PCR
    adata.obs['random'] = np.random.normal(size=adata.n_obs)

    try:
        return scib.metrics.pcr(adata, covariate='random')
    except ZeroDivisionError:
        return float('nan')
    

def pcr_batch(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    try:
        return scib.metrics.pcr(adata, covariate=batch_key)
    except ZeroDivisionError:
        return float('nan')
    

def pcr_label(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    try:
        return scib.metrics.pcr(adata, covariate=label_key)
    except ZeroDivisionError:
        return float('nan')