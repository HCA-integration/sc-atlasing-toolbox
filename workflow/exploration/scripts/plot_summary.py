import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


dfs = [pd.read_table(x) for x in snakemake.input.tsv]
df = pd.concat(dfs).set_index('dataset')
df.to_csv(snakemake.output.tsv, sep='\t')

# plot per dataset
n_rows = df.shape[0]
df.sort_values(by='n_cells', ascending=True).plot.barh(
    subplots=True,
    sharex=False,
    sharey=True,
    layout=(int(df.shape[1] / 2), 2),
    figsize=(12, np.max([5, 1.2 * n_rows]))
)

plt.tight_layout()
plt.savefig(snakemake.output.png)

# save aggregated stats
df[['n_cells', 'n_samples', 'n_donors']].agg('sum').to_csv(snakemake.output.aggregate, sep='\t')