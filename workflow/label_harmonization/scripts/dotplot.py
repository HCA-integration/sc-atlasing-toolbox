import logging
logging.basicConfig(level=logging.INFO)
from typing import MutableMapping
from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense


def filter_markers(markers, _keys):
    if isinstance(markers, dict):
        markers = {
            k: [g for g in _keys if g in v]
            for k, v in markers.items()
        }
        return {k: l for k, l in markers.items() if len(l) > 0}
    elif isinstance(markers, list):
        return [g for g in markers if g in _keys]
    else:
        raise ValueError('marker must be dict or list')


def mapping_to_dict(mapping: MutableMapping) -> dict:
    return {} if mapping is None else dict(mapping.items())


input_file = snakemake.input[0]
output_file = snakemake.output[0]
group = snakemake.wildcards.group

kwargs = mapping_to_dict(snakemake.params.kwargs)
marker_genes = snakemake.params.marker_genes

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    obs='obs',
    var='var',
    X='X',
    dask=True,
    backed=True,
)

# match marker genes and var_names
logging.info(adata.var)
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

# filter marker genes
logging.info('Marker genes before filtering:')
logging.info(marker_genes)
marker_genes = filter_markers(marker_genes, adata.var_names)
assert len(marker_genes) > 0, 'No marker genes found in data'
logging.info('Marker genes after filtering:')
logging.info(marker_genes)
n_marker_genes = len(marker_genes) if isinstance(marker_genes, list) else sum(len(v) for v in marker_genes.values())
logging.info(f'Number of marker genes: {n_marker_genes}')

# subset to genes for plotting
if isinstance(marker_genes, dict):
    genes_to_plot = set([element for sublist in marker_genes.values() for element in sublist])
    genes_to_plot = list(genes_to_plot)
else:
    genes_to_plot = marker_genes
logging.info(f'Load {len(genes_to_plot)} genes to memory...')
adata = adata[:, adata.var_names.isin(genes_to_plot)].copy()
adata = dask_compute(adata)
adata = ensure_dense(adata)
logging.info(adata.__str__())

if group == 'all':
    logging.info('Plotting global dotplot...')
    sc.pl.dotplot(
        adata,
        groupby='group',
        var_names=marker_genes,
        show=False,
        title=f'Markers vs groups for {snakemake.wildcards.dataset}',
        swap_axes=False,
        # swap_axes=n_marker_genes > adata.obs['group'].nunique(),
        **kwargs,
    )
    plt.savefig(output_file, bbox_inches='tight', dpi=200)
    exit()

adata = adata[adata.obs['group'] == group]
if adata.n_obs == 0:
    logging.info('No cells, skip...')
    plt.savefig(output_file, bbox_inches='tight', dpi=200)
    exit()

logging.info(f'Plotting for group={group}...')
sc.pl.dotplot(
    adata,
    groupby=['reannotation_index', 'reannotation'],
    var_names=filter_markers(marker_genes, adata.var_names),
    show=False,
    title=f'Group: {group}',
    **kwargs,
)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
