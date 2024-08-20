import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import cellhint

dpi = 100
plt.rcParams['figure.dpi'] = dpi

input_file = snakemake.input[0]
output_treeplot = Path(snakemake.output.treeplot)
output_treeplot_ordered = Path(snakemake.output.treeplot_ordered)
output_heatmap = snakemake.output.heatmap
# output_sankeyplot = snakemake.output.sankeyplot
dataset = snakemake.wildcards.dataset

output_treeplot.mkdir(parents=True, exist_ok=True)
output_treeplot_ordered.mkdir(parents=True, exist_ok=True)

logging.info(f'Load cell type alignment from {input_file}...')
alignment = cellhint.DistanceAlignment.load(input_file)

groups = list(np.unique(alignment.groups))
for group in groups + ['all']:
    logging.info(f'Treeplots or group={group}...')
    try:
        if group == 'all':
            algn = alignment.relation
        else:
            algn = alignment.relation[alignment.groups == group]
            print(algn, flush=True)
        cellhint.treeplot(
            algn,
            order_dataset=False,
            label_size=12,
            # figsize = (20, 30),
            title=f'CellTypist label harmonization: {dataset}',
            show=False,
        )
        plt.savefig(output_treeplot / f'{group}.png', dpi=dpi, bbox_inches='tight')

        cellhint.treeplot(
            algn,
            order_dataset=True,
            label_size=12,
            # figsize = (20, 30),
            title=f'CellTypist label harmonization: {dataset}',
            show=False,
            # save=output_treeplot_ordered / f'{group}.png'
        )
        plt.savefig(output_treeplot_ordered / f'{group}.png', dpi=dpi, bbox_inches='tight')
    except Exception:
        logging.error(f'Saving empty treeplots for group={group}')
        plt.figure(figsize=(6, 4))
        plt.savefig(output_treeplot / f'{group}.png')
        plt.savefig(output_treeplot_ordered / f'{group}.png')

logging.info('Heatmap...')
plot = sns.clustermap(alignment.base_distance.to_meta())
plot.fig.suptitle(dataset)
plt.savefig(output_heatmap, dpi=dpi, bbox_inches='tight')

# cellhint.sankeyplot(alignment, show=False, save=output_sankeyplot)