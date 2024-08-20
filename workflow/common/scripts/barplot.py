import numpy as np
import pandas as pd
import seaborn as sns

input_tsv = snakemake.input.tsv
output_png = snakemake.output.png


params = snakemake.params
category = params.get('category')
metric = params.get('metric')
facet_row = params.get('facet_row')
facet_col = params.get('facet_col')
hue = params.get('hue')

xlim = params.get('xlim')
ylim = params.get('ylim')

title = params.get("title", f"{metric} for {category}")
description = params.get("description", "")

# read table
df = pd.read_table(input_tsv)

# convert to category
df[category] = df[category].astype('category')
for col in [facet_row, facet_col, hue]:
    if col is not None:
        df[col] = df[col].astype(str).astype('category')


if hue in df.columns:
    hue = None if df[hue].nunique() > 6 else hue

# plot parameters
n_rows = 1 if facet_row is None else df[facet_row].nunique()
n_cols = 1 if facet_col is None else df[facet_col].nunique()
adjust = np.min([.8 + (.04 * n_rows), .95])
order = df.groupby(category, observed=False)[metric].min().sort_values(ascending=False).index

g = sns.catplot(
    data=df,
    x=metric,
    y=category,
    row=facet_row,
    col=facet_col,
    hue=hue,
    palette='Paired' if hue else None,
    margin_titles=True,
    kind='bar',
    order=order,
    err_kws={'linewidth': 0.5},
    capsize=0.2,
    dodge=params.get('dodge', False),
    aspect=1.5,
    legend_out=False,
)
g.set(xlim=xlim, ylim=ylim)

# configure legend
if hue:
    sns.move_legend(g, 'upper center', bbox_to_anchor=(0.6, 0))

g.fig.subplots_adjust(top=adjust)
g.fig.suptitle(f'{title} {description}')

g.savefig(output_png, dpi=200, bbox_inches='tight')
