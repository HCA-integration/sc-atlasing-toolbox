from matplotlib import pyplot as plt
import seaborn as sns
import celltypist

dpi = 300
plt.rcParams['figure.dpi'] = dpi

input_file = snakemake.input[0]
output_treeplot = snakemake.output.treeplot
output_treeplot_ordered = snakemake.output.treeplot_ordered
output_heatmap = snakemake.output.heatmap
# output_sankeyplot = snakemake.output.sankeyplot
dataset = snakemake.wildcards.dataset

alignment = celltypist.DistanceAlignment.load(input_file)

celltypist.treeplot(
    alignment,
    order_dataset=False,
    label_size=12,
    # figsize = (20, 30),
    title=f'CellTypist label harmonization: {dataset}',
    show=False,
    save=output_treeplot
)

celltypist.treeplot(
    alignment,
    order_dataset=True,
    label_size=12,
    # figsize = (20, 30),
    title=f'CellTypist label harmonization: {dataset}',
    show=False,
    save=output_treeplot_ordered
)

plot = sns.clustermap(alignment.base_distance.to_meta())
plot.fig.suptitle(dataset)
plt.savefig(output_heatmap, dpi=dpi, bbox_inches='tight')

# celltypist.sankeyplot(alignment, show=False, save=output_sankeyplot)