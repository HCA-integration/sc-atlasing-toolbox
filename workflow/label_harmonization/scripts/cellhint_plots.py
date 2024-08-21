import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import cellhint


input_file = snakemake.input[0]
output_treeplot = snakemake.output.treeplot
output_treeplot_ordered = snakemake.output.treeplot_ordered
output_heatmap = snakemake.output.heatmap
# output_sankeyplot = snakemake.output.sankeyplot

dataset = snakemake.wildcards.dataset
group = snakemake.wildcards.group


logging.info(f'Load cell type alignment from {input_file}...')
alignment = cellhint.DistanceAlignment.load(input_file)

# subset group
if group == 'all':
    relation = alignment.relation
    dist_mat = alignment.base_distance.to_meta()
else:
    relation = alignment.relation[alignment.groups == group]
    print('relation:\n', relation, flush=True)
    sub_reannotations = alignment.reannotation.query('group == @group')[['dataset', 'cell_type']].drop_duplicates().agg(': '.join, axis=1)
    print('reannotations:\n', sub_reannotations, flush=True)
    dist_mat = alignment.base_distance.to_meta().loc[sub_reannotations, sub_reannotations]
    print('dist_mat:\n', dist_mat, flush=True)

dpi = 100
plt.rcParams['figure.dpi'] = dpi

logging.info('Heatmap...')
if dist_mat.shape[0] > 1:
    g = sns.clustermap(dist_mat)
    g.fig.suptitle(dataset)
    plt.savefig(output_heatmap, dpi=dpi, bbox_inches='tight')
else:
    logging.error(f'Empty distance matrix for group={group}')
    plt.figure(figsize=(6, 4))
    plt.savefig(output_heatmap)

logging.info(f'Treeplots or group={group}...')
try:
    cellhint.treeplot(
        relation,
        order_dataset=False,
        label_size=12,
        # figsize = (20, 30),
        title=f'CellHint label harmonization: {dataset}',
        show=False,
    )
    plt.savefig(output_treeplot, dpi=dpi, bbox_inches='tight')

    cellhint.treeplot(
        relation,
        order_dataset=True,
        label_size=12,
        # figsize = (20, 30),
        title=f'CellHint label harmonization: {dataset}',
        show=False,
        # save=output_treeplot_ordered
    )
    plt.savefig(output_treeplot_ordered, dpi=dpi, bbox_inches='tight')
except Exception:
    logging.error(f'Saving empty treeplots for group={group}')
    plt.figure(figsize=(6, 4))
    plt.savefig(output_treeplot)
    plt.savefig(output_treeplot_ordered)


# cellhint.sankeyplot(alignment, show=False, save=output_sankeyplot)