import anndata as ad
import scanpy as sc
import numpy as np

from utils.accessors import subset_hvg

def morans_i_gene_score(_ad, _gene_signature, var_key) -> float:
    # Subset to highly variable genes
    _ad.var[var_key] = _ad.var[var_key]
    _ad, _ = subset_hvg(_ad, var_column=var_key, min_cells=1, compute_dask=True)


    # TODO: how to subset for hvg AND signature?
    # --> currently only using hvg 
    
    #_ad.var['subset_genes'] = _ad.var['highly_variable'] | _ad.var_names.isin(ifn_signature)
    #_ad.var['subset_genes'].value_counts()
    _ad.var[var_key] = _ad.var[var_key]
    _ad, _ = subset_hvg(_ad, var_column=var_key, min_cells=1, compute_dask=True)

    _gene_signature = [g for g in _gene_signature if g in _ad.var_names]

    # Create gene score for the gene signature
    sc.tl.score_genes(_ad, gene_list=_gene_signature, score_name='gene_score')

    try:
        return sc.metrics.morans_i(_ad, vals=_ad.obs['gene_score'])
    except ZeroDivisionError:
        return float('nan')

# Moran's I for platlets -> baseline
def morans_i_platlets(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    return morans_i_gene_score(adata, ["GP1BB", "ITGA2B", "PF4", "PPBP", "TUBB1"], var_key)

# Moran's I for red blood cells -> baseline
def morans_i_rbc(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    return morans_i_gene_score(adata, ["HBA1", "HBA2", "HBB"], var_key)

# Moran's I for plasma cells
def morans_i_plasma_cells(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    return morans_i_gene_score(adata, ["CD79A", "DERL3", "IGHA1", "ITM2C", "JCHAIN", "MZB1", "POU2AF1", "TNFRSF17", "TXNDC11"], var_key)

# Moran's I for T CD4 CTL cells
def morans_i_t_cd4_ctl(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    adata = ad.AnnData(X=adata_raw.X, var=adata_raw.var, obs=adata.obs, obsp=adata.obsp, uns=adata.uns)

    # Yoshida specific?
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    return morans_i_gene_score(adata, ["CCL5", "FGFBP2", "GZMA", "GNLY", "GZMB", "GZMH", "ITGB1", "KLRB1"], var_key)