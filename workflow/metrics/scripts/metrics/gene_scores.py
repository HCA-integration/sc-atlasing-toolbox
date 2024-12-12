import anndata as ad
import scanpy as sc
import numpy as np
import scib

from utils.accessors import subset_hvg

dict_marker_genes = {
    "ifn_signature": ["IRF7", "XAF1", "UBE2L6", "TRIM22", "STAT1",
                        "SP110", "SAMD9L", "SAMD9", "PLSCR1", "PARP9",
                        "OAS2", "OAS1", "MX2", "MX1", "LY6E",
                        "ISG15", "IFIT3", "IFI6", "IFI44L", "IFI35",
                         "HERC5", "EPSTI1", "EIF2AK2", "CMPK2", "BST2"],
    "platlets": ["GP1BB", "ITGA2B", "PF4", "PPBP", "TUBB1"],
    "rbc": ["HBA1", "HBA2", "HBB"],
    "plasma_cells": ["CD79A", "DERL3", "IGHA1", "ITM2C", "JCHAIN", "MZB1", "POU2AF1", "TNFRSF17", "TXNDC11"],
    "t_cd4_ctl": ["CCL5", "FGFBP2", "GZMA", "GNLY", "GZMB", "GZMH", "ITGB1", "KLRB1"],
    "nk_cells": ["GNLY", "KLRB1", "KLRD1", "NCAM1", "NCR1", "NKG7", "TYROBP"],
    "ilc":["KIT", "KLRB1", "LINC01229", "TNFRSF4", "TNFRSF18", "TRDC", "TTLL10"],
    "monocytes_cd14":["CD14", "CSF3R", "S100A12", "S100A8", "S100A9"]
}

def get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key) -> float:
    adata = ad.AnnData(
        X=adata_raw.X,
        var=adata_raw.var,
        obs=adata.obs,
        obsp=adata.obsp,
        uns=adata.uns
    )
    adata.var_names = adata.var_names.astype(str)

    gene_signature = [g for g in gene_signature if g in adata.var_names]
    adata.var.loc[adata.var_names.isin(gene_signature), var_key] = True

    adata, _ = subset_hvg(adata, var_column=var_key, min_cells=1, compute_dask=True)
    
    sc.tl.score_genes(adata, gene_list=gene_signature, score_name='gene_score')
    return adata


def get_bootstrap_adata(adata, size):
    rand_indices = np.random.choice(adata.obs.index, size=size, replace=False)
    adata_bootstrap = adata[adata.obs.index.isin(rand_indices)]
    # Recalculation is important as kNN graph is input for Moran's I
    try: 
        adata_bootstrap.obsm['distances'] = adata_bootstrap.obsp['distances']
        sc.pp.neighbors(adata_bootstrap, use_rep='distances', metric='precomputed', transformer='sklearn')   
    except ValueError:
        sc.pp.neighbors(adata_bootstrap, use_rep='X_emb')   
    finally:  
        return adata_bootstrap


def get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, size):
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)

    morans_values = []
    for _ in range(5):
        ad = get_bootstrap_adata(adata, size)
        morans_values.append(sc.metrics.morans_i(ad, vals=ad.obs['gene_score']))
    
    return morans_values


# Moran's I for platlets -> baseline
def morans_i_platlets(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("platlets")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for red blood cells -> baseline
def morans_i_rbc(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("rbc")

    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for plasma cells
def morans_i_plasma_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("plasma_cells")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for T CD4 CTL cells
def morans_i_t_cd4_ctl(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("t_cd4_ctl")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for NK cells
def morans_i_nk_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("nk_cells")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for ILC cells
def morans_i_ilc_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("ilc")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran's I for Monocyte CD14 cells
def morans_i_monocytes_cd14(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("monocytes_cd14")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Moran' I for the IFN gene signature (specific for the data set by Yoshida et al.)
def morans_i_ifn_signature(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("ifn_signature")
    
    return get_morans_i_gene_score_bt(adata, adata_raw, gene_signature, var_key, int(0.9*adata.n_obs))


# Principle Component Regression for platlets -> baseline
def pcr_platlets(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("platlets")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for red blood cells -> baseline
def pcr_rbc(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("rbc")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for plasma cells
def pcr_plasma_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("plasma_cells")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for T CD4 CTL cells
def pcr_t_cd4_ctl(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("t_cd4_ctl")

    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for NK cells
def pcr_nk_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("nk_cells")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for ILC cells
def pcr_ilc_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("ilc")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for Monocyte CD14 cells
def pcr_monocytes_cd14(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("monocytes_cd14")
    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)


# Principle Component Regression for the IFN gene signature (specific for the data set by Yoshida et al.)
def pcr_ifn_signature(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    gene_signature = dict_marker_genes.get("ifn_signature")

    adata = get_anndata_for_gene_score(adata, adata_raw, gene_signature, var_key)
    
    return scib.metrics.pcr(adata, covariate='gene_score', recompute_pca=False)