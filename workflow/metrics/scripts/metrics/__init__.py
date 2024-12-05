from .ari import ari, ari_kmeans_y, ari_leiden_y
from .asw import asw_batch, asw_batch_y, asw_label, asw_label_y
from .graph_connectivity import graph_connectivity, graph_connectivity_y
from .isolated_labels import isolated_label_asw, isolated_label_f1, isolated_label_asw_y
from .lisi import clisi, clisi_y, ilisi, ilisi_y
from .nmi import nmi, nmi_kmeans_y, nmi_leiden_y
from .pcr import cell_cycle, pcr, pcr_y, pcr_batch, pcr_label, pcr_random
from .kbet import kbet_y
from .morans_i import morans_i_batch, morans_i_label, morans_i_random
from .gene_scores import *

metric_map = {
    'ari': ari,
    'ari_kmeans_y': ari_kmeans_y,
    'ari_leiden_y': ari_leiden_y,
    'asw_batch': asw_batch,
    'asw_batch_y': asw_batch_y,
    'asw_label': asw_label,
    'asw_label_y': asw_label_y,
    'asw_label_y': asw_label_y,
    'graph_connectivity': graph_connectivity,
    'graph_connectivity_y': graph_connectivity_y,
    'isolated_label_asw': isolated_label_asw,
    'isolated_label_asw_y': isolated_label_asw_y,
    'isolated_label_f1': isolated_label_f1,
    'clisi': clisi,
    'clisi_y': clisi_y,
    'ilisi': ilisi,
    'ilisi_y': ilisi_y,
    'clisi': clisi,
    'clisi': clisi,
    'nmi': nmi,
    'nmi_kmeans_y': nmi_kmeans_y,
    'nmi_leiden_y': nmi_leiden_y,
    'cell_cycle': cell_cycle,
    'pcr': pcr,
    'pcr_y': pcr_y,
    'kbet_y': kbet_y,
    'morans_i_random': morans_i_random,
    'morans_i_batch': morans_i_batch,
    'morans_i_label': morans_i_label,
    'morans_i_platlets': morans_i_platlets,
    'morans_i_rbc': morans_i_rbc,
    'morans_i_plasma_cells': morans_i_plasma_cells,
    'morans_i_t_cd4_ctl': morans_i_t_cd4_ctl,
    'morans_i_nk_cells': morans_i_nk_cells,
    'morans_i_ifn_signature': morans_i_ifn_signature,
    'pcr_random': pcr_random,
    'pcr_batch': pcr_batch,
    'pcr_label': pcr_label,
    'pcr_platlets': pcr_platlets,
    'pcr_rbc': pcr_rbc,
    'pcr_plasma_cells': pcr_plasma_cells,
    'pcr_t_cd4_ctl': pcr_t_cd4_ctl,
    'pcr_nk_cells': pcr_nk_cells,
    'pcr_ifn_signature': pcr_ifn_signature,
}