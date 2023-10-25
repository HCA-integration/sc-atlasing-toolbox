import numpy as np
import pandas as pd
import seaborn as sns


input_tsv = snakemake.input.tsv
output_png = snakemake.output.png

params_keys = [
    'metric',
    'category',
    'hue',
    'facet_row',
    'facet_col',
    'title',
    'description',
    'xlim',
    'ylim'
]
params = {p: None for p in params_keys}
params.update(snakemake.params)

category = params['category']
metric = params['metric']
facet_row = params['facet_row']
facet_col = params['facet_col']

df = pd.read_table(input_tsv)

# convert to category
df[category] = df[category].astype('category')
for col in [facet_row, facet_col]:
    if col is not None:
        df[col] = df[col].astype('category')

if hue in df.columns:
    hue = None if df[hue].nunique() > 6 else hue

# plot parameters
n_rows = 1 if facet_row is None else df[facet_row].nunique()
n_cols = 1 if facet_col is None else df[facet_col].nunique()
n_cats = df[category].nunique()
adjust = np.min([.8 + (.04 * n_rows), .95])
order = df.groupby(category)[metric].min().sort_values(ascending=False).index
aspect = np.max([1.5, (.15 * n_cats)])

g = sns.catplot(
    data=df,
    x=category,
    y=metric,
    sharey=False,
    row=facet_row,
    col=facet_col,
    hue=params['hue'],
    height=3,
    aspect=aspect,
    margin_titles=True,
    kind='swarm',
    s=10,
    order=order,
)
g.set(ylim=params['ylim'])
g.tick_params(axis='x', rotation=90)
g.fig.subplots_adjust(top=adjust)
g.fig.suptitle(f'{params["title"]} {params["description"]}')

g.savefig(output_png)
