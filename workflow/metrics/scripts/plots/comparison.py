import numpy as np
import pandas as pd
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')


metrics_df = snakemake.input.metrics
metrics_df['metric'] = metrics_df['metric'].astype('category')

theislab_metrics = ['nmi', 'ari', 'ilisi', 'clisi', 'asw_label', 'asw_batch']
yoseflab_metrics = ['nmi_leiden_y', 'ari_leiden_y', 'ilisi_y', 'clisi_y', 'asw_label_y', 'asw_batch_y']

metrics_df.loc[metrics_df['metric'].isin(theislab_metrics), 'lab'] = 'theislab'
metrics_df.loc[metrics_df['metric'].isin(yoseflab_metrics), 'lab'] = 'yoseflab'

metric_name_map = {
    'ari_leiden_y': 'ari',
    'nmi_leiden_y': 'nmi',
    'ilisi_y': 'ilisi',
    'clisi_y': 'clisi',
    'asw_label_y': 'asw_label',
    'asw_batch_y': 'asw_batch'
}

metrics_df['metric'] = [metric_name_map[x] if x in metric_name_map.keys() else x for x in metrics_df['metric']]

df = metrics_df.pivot(
    index=['dataset', 'batch', 'label', 'hyperparams', 'metric', 'method', 'output_type'],
    columns='lab',
    values='s'
).reset_index()
df['speedup'] = df['theislab'] / df['yoseflab']

g = sns.catplot(
    data=df,
    x='speedup',
    y='metric',
    #row=facet_var,
    col='output_type',
    hue='output_type',
    #height=3,
    #aspect=0.8 * n_cols,
    margin_titles=True,
    kind='bar',
    order=metrics_df.groupby('metric')['s'].max().sort_values(ascending=False).index,
    errwidth=0.5,
    capsize=0.2,
    dodge=False,
)
