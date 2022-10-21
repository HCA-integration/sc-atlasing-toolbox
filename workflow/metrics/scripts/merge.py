import numpy as np
import pandas as pd
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

wildcards = snakemake.wildcards
expanded_wildcards = snakemake.params['wildcards']
expanded_wildcards = pd.DataFrame(expanded_wildcards)
facet_var = expanded_wildcards.columns[0]

metrics_df = pd.concat([pd.read_table(file) for file in snakemake.input.metrics]).reset_index(drop=True)
benchmark_df = pd.concat([pd.read_table(file) for file in snakemake.input.benchmark]).reset_index(drop=True)

benchmark_df = pd.concat([expanded_wildcards, benchmark_df], axis=1)

metrics_df = metrics_df.merge(benchmark_df)
metrics_df.to_csv(snakemake.output.metrics, sep='\t', index=False)

print(metrics_df[['metric', 'method', 'output_type', 'metric_type', 'score', 's', 'h:m:s']])

## Plots
metrics_df['metric'] = metrics_df['metric'].astype('category')
metrics_df[facet_var] = metrics_df[facet_var].astype('category')
n_cols = metrics_df['metric_type'].nunique()
n_rows = metrics_df[facet_var].nunique()

# metrics plot
g = sns.catplot(
    data=metrics_df,
    x='metric',
    y='score',
    row=facet_var,
    col='output_type',
    hue='metric_type',
    height=3,
    aspect=0.8 * n_cols,
    margin_titles=True,
    s=10,
    order=metrics_df.groupby('metric')['score'].max().sort_values(ascending=False).index,
)
g.tick_params(axis='x', rotation=90)
description = ' '.join([f'{key}={value}' for key, value in wildcards.items()])
g.fig.subplots_adjust(top=.8 + (.04 * n_rows))
g.fig.suptitle(f'Metric scores {description}')
g.savefig(snakemake.output.plot)


# Time plot
g = sns.catplot(
    data=metrics_df,
    x='s',
    y='metric',
    row=facet_var,
    col='output_type',
    hue='metric_type',
    height=3,
    aspect=0.8 * n_cols,
    margin_titles=True,
    kind='bar',
    order=metrics_df.groupby('metric')['s'].max().sort_values(ascending=False).index,
    errwidth=0.5,
    capsize=0.2,
    dodge=False,
)
description = ' '.join([f'{key}={value}' for key, value in wildcards.items()])
g.fig.subplots_adjust(top=.8 + (.04 * n_rows))
g.fig.suptitle(f'Computation time for metrics {description}')
g.savefig(snakemake.output.time)
