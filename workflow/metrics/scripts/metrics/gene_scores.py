import anndata as ad
import scanpy as sc
import numpy as np
import scib

from utils.accessors import subset_hvg
from .bootstrap import bootstrap_metric
from .morans_i import _morans_i


# TODO: get from user config file
dict_marker_genes = {
    # IFN gene signature (specific for the data set by Yoshida et al.)
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


def get_gene_score(adata, adata_raw, gene_signature, var_key):
    adata = ad.AnnData(
        X=adata_raw.X,
        var=adata_raw.var,
        obs=adata.obs,
        obsp=adata.obsp,
        uns=adata.uns
    )

    gene_signature = [g for g in gene_signature if g in adata.var_names]
    # also include genes of interest to subset
    adata.var.loc[adata.var_names.isin(gene_signature), var_key] = True

    # subset to HVGs (used for control genes) + genes of interest
    adata, _ = subset_hvg(adata, var_column=var_key, min_cells=1, compute_dask=True)
    
    # calculate gene score
    sc.tl.score_genes(adata, gene_list=gene_signature)
    
    return adata.obs['score']


def morans_i_gene_score(
    adata,
    adata_raw,
    var_key,
    gene_signature,
    n_bootstraps=5,
    bootstrap_size=None,
):
    """
    Moran's I for gene score
    """
    try :
        gene_scores = get_gene_score(adata, adata_raw, gene_signature, var_key)
    except ValueError:
        return float('nan')
        
    return bootstrap_metric(
        adata,
        metric_function=lambda _ad, vals: _morans_i(_ad, vals=vals[_ad.obs_names]),
        n_bootstraps=n_bootstraps,
        size=bootstrap_size,
        vals=gene_scores,
    )


def morans_i_random_genes(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    random_genes = np.random.choice(adata_raw.var_names[adata_raw.var[var_key]], size=50, replace=False)
    return morans_i_gene_score(adata, adata_raw, var_key, random_genes)


def morans_i_platlets(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["platlets"])


# Moran's I for red blood cells -> baseline
def morans_i_rbc(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["rbc"])


# Moran's I for plasma cells
def morans_i_plasma_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["plasma_cells"])


# Moran's I for T CD4 CTL cells
def morans_i_t_cd4_ctl(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["t_cd4_ctl"])


# Moran's I for NK cells
def morans_i_nk_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["nk_cells"])


def morans_i_ilc_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["ilc"])


def morans_i_monocytes_cd14(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["monocytes_cd14"])


def morans_i_ifn_signature(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return morans_i_gene_score(adata, adata_raw, var_key, dict_marker_genes["ifn_signature"])


def pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature):
    """   
    Principle Component Regression for gene score
    """
    if output_type == 'knn':
        return float('nan')
    
    try:
        adata.obs['gene_score'] = get_gene_score(adata, adata_raw, gene_signature, var_key)
    except ValueError:
        return float('nan')
    
    return scib.metrics.pcr(
        adata,
        covariate='gene_score',
        recompute_pca=False
    )


def pcr_random_genes(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    random_genes = np.random.choice(adata_raw.var_names[adata_raw.var[var_key]], size=50, replace=False)
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=random_genes)


def pcr_platlets(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['platlets'])


def pcr_rbc(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['rbc'])


def pcr_plasma_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['plasma_cells'])


def pcr_t_cd4_ctl(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['t_cd4_ctl'])


def pcr_nk_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['nk_cells'])


def pcr_ilc_cells(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['ilc'])


def pcr_monocytes_cd14(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['monocytes_cd14'])


def pcr_ifn_signature(adata, output_type, adata_raw, var_key='metrics_features', **kwargs):
    return pcr_gene_score(adata, output_type, adata_raw, var_key, gene_signature=dict_marker_genes['ifn_signature'])