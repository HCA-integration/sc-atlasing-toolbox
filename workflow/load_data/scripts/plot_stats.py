import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


in_stats = snakemake.input
out_intersection_stats = snakemake.output.intersection_stats
out_intersection_plot = snakemake.output.intersection_plot

series = [
    pd.read_table(f, index_col=0, header=0, names=[study]).T
    for study, f in in_stats.items()
]
stats_df = pd.concat(series, axis=0)
stats_df.index.name = 'study'
stats_df['intersection_fraction'] = stats_df['intersection_fraction'].astype(float).round(2)
stats_df.reset_index().to_csv(out_intersection_stats, sep='\t', index=False)

plot = stats_df['intersection_fraction'].sort_values().plot.barh(
    title='Intersection of donor_id with best matching DCP ID'
)
for container in plot.containers:
    plot.bar_label(container)
plot.set_xlabel('ID overlap fraction (normed over CxG donors/samples)')
plt.savefig(out_intersection_plot, bbox_inches='tight')