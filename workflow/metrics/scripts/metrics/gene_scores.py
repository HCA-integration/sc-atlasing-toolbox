import anndata as ad
import scanpy as sc
import numpy as np
import scib

from utils.accessors import subset_hvg

def get_anndata_for_gene_score(_adata, _adata_raw, _gene_signature, _var_key) -> float:
    _ad = ad.AnnData(X=_adata_raw.X, var=_adata_raw.var, obs=_adata.obs, obsp=_adata.obsp, uns=_adata.uns)

    # Subset to highly variable genes and the gene signature
    _gene_signature = [g for g in _gene_signature if g in _ad.var_names]
    _ad.var[_var_key] = _ad.var[_var_key]

    _ad.var.loc[_ad.var_names.isin(_gene_signature), _var_key] = True
    _ad, _ = subset_hvg(_ad, var_column=_var_key, min_cells=1, compute_dask=True)

    # Create gene score for the gene signature (only if gene signature is not empty)
    if not _gene_signature:
        return None
    else:
        sc.tl.score_genes(_ad, gene_list=_gene_signature, score_name='gene_score')
        return _ad

# Moran's I for platlets -> baseline
def morans_i_platlets(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["GP1BB", "ITGA2B", "PF4", "PPBP", "TUBB1"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)

    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score']) # if adata is not None else float('nan')


# Moran's I for red blood cells -> baseline
def morans_i_rbc(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["HBA1", "HBA2", "HBB"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score'])


# Moran's I for plasma cells
def morans_i_plasma_cells(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["CD79A", "DERL3", "IGHA1", "ITM2C", "JCHAIN",
                      "MZB1", "POU2AF1", "TNFRSF17", "TXNDC11"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score'])


# Moran's I for T CD4 CTL cells
def morans_i_t_cd4_ctl(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["CCL5", "FGFBP2", "GZMA", "GNLY", "GZMB",
                      "GZMH", "ITGB1", "KLRB1"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score'])


# Moran's I for NK cells
def morans_i_nk_cells(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["GNLY", "KLRB1", "KLRD1", "NCAM1", "NCR1",
                      "NKG7", "TYROBP"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score'])


# Moran' I for the IFN gene signature (specific for the data set by Yoshida et al.)
def morans_i_ifn_signature(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["IRF7", "XAF1", "UBE2L6", "TRIM22", "STAT1",
                      "SP110", "SAMD9L", "SAMD9", "PLSCR1", "PARP9",
                      "OAS2", "OAS1", "MX2", "MX1", "LY6E",
                      "ISG15", "IFIT3", "IFI6", "IFI44L", "IFI35",
                      "HERC5", "EPSTI1", "EIF2AK2", "CMPK2", "BST2"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return sc.metrics.morans_i(adata, vals=adata.obs['gene_score'])


# Principle Component Regression for platlets -> baseline
def pcr_platlets(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["GP1BB", "ITGA2B", "PF4", "PPBP", "TUBB1"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for red blood cells -> baseline
def pcr_rbc(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["HBA1", "HBA2", "HBB"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for plasma cells
def pcr_plasma_cells(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["CD79A", "DERL3", "IGHA1", "ITM2C", "JCHAIN",
                      "MZB1", "POU2AF1", "TNFRSF17", "TXNDC11"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for T CD4 CTL cells
def pcr_t_cd4_ctl(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["CCL5", "FGFBP2", "GZMA", "GNLY", "GZMB",
                      "GZMH", "ITGB1", "KLRB1"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for NK cells
def pcr_nk_cells(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["GNLY", "KLRB1", "KLRD1", "NCAM1", "NCR1",
                      "NKG7", "TYROBP"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for the IFN gene signature (specific for the data set by Yoshida et al.)
def pcr_ifn_signature(adata, output_type, batch_key, label_key, adata_raw, var_key='metrics_features', n_threads=1, **kwargs):
    gene_signature = ["IRF7", "XAF1", "UBE2L6", "TRIM22", "STAT1",
                      "SP110", "SAMD9L", "SAMD9", "PLSCR1", "PARP9",
                      "OAS2", "OAS1", "MX2", "MX1", "LY6E",
                      "ISG15", "IFIT3", "IFI6", "IFI44L", "IFI35",
                      "HERC5", "EPSTI1", "EIF2AK2", "CMPK2", "BST2"]
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)