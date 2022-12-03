import numpy as np
import pandas as pd
import seaborn as sns


input_tsv = snakemake.input.tsv
output_png = snakemake.output.png
metric = snakemake.params.metric
category = snakemake.params.category
hue = snakemake.params.hue
facet_row = snakemake.params.facet_row
facet_col = snakemake.params.facet_col
dodge = snakemake.params.dodge
title = snakemake.params.title
description = snakemake.params.description

df = pd.read_table(input_tsv)

# convert to category
df[category] = df[category].astype('category')
for col in [facet_row, facet_col]:
    if col is not None:
        df[col] = df[col].astype('category')

# plot parameters
n_rows = 1 if facet_row is None else df[facet_row].nunique()
n_cols = 1 if facet_col is None else df[facet_col].nunique()
adjust = np.min([.8 + (.04 * n_rows), .95])
order = df.groupby(category)[metric].min().sort_values(ascending=False).index

g = sns.catplot(
    data=df,
    x=metric,
    y=category,
    row=facet_row,
    col=facet_col,
    hue=hue,
    margin_titles=True,
    kind='bar',
    order=order,
    errwidth=0.5,
    capsize=0.2,
    dodge=dodge,
)
g.set(xlim=(-.01, None))
g.fig.subplots_adjust(top=adjust)
g.fig.suptitle(f'{title} {description}')

g.savefig(output_png)
