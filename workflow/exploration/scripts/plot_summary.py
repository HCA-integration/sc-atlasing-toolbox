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
    layout=(2, 2),
    figsize=(12, 2 * n_rows)
)
plt.tight_layout()
plt.savefig(snakemake.output.png)
