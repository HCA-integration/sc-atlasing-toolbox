from matplotlib import pyplot as plt
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense

input_file = snakemake.input[0]
output_rankplot = snakemake.output.rankplot
output_dotplot = snakemake.output.dotplot
output_heatmap = snakemake.output.heatmap

wildcards = snakemake.wildcards
group_key = wildcards.group
file_id = wildcards.file_id
title = f'Marker genes for file_id={file_id} group={group_key}'

args = snakemake.params.get('args', {})
args['key'] = f'marker_genes_group={group_key}'

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

# subset to genes for plotting
marker_gene_key = args.get('key')
genes_to_plot = list(set().union(*adata.uns[marker_gene_key]['names']))
logging.info(f'Load {len(genes_to_plot)} genes to memory...')
adata = adata[:, adata.var_names.isin(genes_to_plot)].copy()
adata = dask_compute(adata)
adata = ensure_dense(adata)
logging.info(adata.__str__())

logging.info('Rankplot...')
sc.pl.rank_genes_groups(adata, **args)
plt.suptitle(title)
plt.savefig(output_rankplot, dpi=100, bbox_inches='tight')

logging.info('Dotplot...')
sc.pl.rank_genes_groups_dotplot(
    adata,
    **args,
    standard_scale='var',
    swap_axes=False,
    title=title,
)
plt.savefig(output_dotplot, dpi=100, bbox_inches='tight')

logging.info('Heatmap...')
sc.pl.rank_genes_groups_heatmap(
    adata,
    **args,
    standard_scale='var',
    swap_axes=True,
    show_gene_labels=len(genes_to_plot) < 300,
)
plt.suptitle(title)
plt.savefig(output_heatmap, dpi=100, bbox_inches='tight')
