import numpy as np
import scanpy as sc
from pathlib import Path
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense

input_file = snakemake.input[0]
output_rankplot = snakemake.output.rankplot
output_dotplot = snakemake.output.dotplot
output_heatmap = snakemake.output.heatmap

Path(output_dotplot).mkdir(parents=True, exist_ok=True)
Path(output_heatmap).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
group_key = wildcards.group
file_id = wildcards.file_id
title = f'Marker genes for file_id={file_id} group={group_key}'

args = snakemake.params.get('args', {})
args['key'] = f'marker_genes_group={group_key}'
n_groups_per_split = args.pop('n_groups_per_split', 10)

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

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']

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

# partition groups to make plots visualisable
groups = adata.uns[f'marker_genes_group={group_key}']['names'].dtype.names
n_splits = max(1, int(len(groups) / n_groups_per_split))
splits = np.array_split(groups, n_splits)

for i, partition_groups in enumerate(splits):
    # partition_groups = list(partition_groups)
    logging.info(f'Dotplot for partition {i}...')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        **args,
        groups=partition_groups,
        standard_scale='var',
        swap_axes=False,
        title=title,
    )
    plt.savefig(
        f'{output_dotplot}/split={i}.png',
        dpi=100,
        bbox_inches='tight',
    )

    logging.info(f'Heatmap for partition {i}...')
    sc.pl.rank_genes_groups_heatmap(
        adata,
        **args,
        groups=partition_groups,
        standard_scale='var',
        swap_axes=True,
        show_gene_labels=len(genes_to_plot) < 300,
    )
    plt.suptitle(title)
    plt.savefig(
        f'{output_heatmap}/split={i}.png',
        dpi=100,
        bbox_inches='tight'
    )
