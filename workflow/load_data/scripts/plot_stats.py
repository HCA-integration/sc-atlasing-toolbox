import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


in_stats = snakemake.input
out_intersection = snakemake.output.intersection

series = [
    pd.read_table(f, index_col=0, header=0, names=[study]).T
    for study, f in in_stats.items()
]
stats_df = pd.concat(series, axis=0)
print(stats_df)

stats_df['intersection_fraction'].plot.barh(
    title='Intersection of donor_id with best matching DCP ID'
)
plt.savefig(out_intersection, bbox_inches='tight')