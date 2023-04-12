from matplotlib import pyplot as plt
import seaborn as sns
import celltypist


input_file = snakemake.input[0]
output_treeplot = snakemake.output.treeplot
output_heatmap = snakemake.output.heatmap
# output_sankeyplot = snakemake.output.sankeyplot
dataset = snakemake.wildcards.dataset

alignment = celltypist.DistanceAlignment.load(input_file)

celltypist.treeplot(
    alignment,
    order_dataset=False,
    label_size=8,
    title=f'CellTypist label harmonization: {dataset}',
    show=False,
    save=output_treeplot
)

plot = sns.clustermap(alignment.base_distance.to_meta())
plot.fig.suptitle(dataset)
plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')

# celltypist.sankeyplot(alignment, show=False, save=output_sankeyplot)