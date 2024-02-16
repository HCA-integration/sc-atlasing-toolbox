from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds


def get_fraction_removed(df, group, key='passed_qc'):
    grouped_counts = df.groupby([group, key], observed=False).size().reset_index(name='Counts')
    grouped_frac = grouped_counts.pivot(index=group, columns=key, values='Counts')
    grouped_frac['fraction_removed'] = (grouped_frac[False] / grouped_frac.sum(axis=1)).round(2)
    return grouped_frac


input_zarr = snakemake.input.zarr
output_plots = snakemake.output.plots

output_plots = Path(output_plots)
output_plots.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')

# If no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('Empty data, skipping plots...')
    exit()

# get parameters
file_id = snakemake.wildcards.file_id
dataset, groups = parse_parameters(adata, snakemake.params, filter_hues=True)
threshold_keys = ['n_counts', 'n_genes', 'percent_mito'] 
thresholds = get_thresholds(
    threshold_keys=threshold_keys,
    autoqc_thresholds=adata.uns['scautoqc_ranges'],
    user_thresholds=snakemake.params.get('thresholds'),
    user_key=file_id,
)
logging.info(f'\n{pformat(thresholds)}')

logging.info('Apply thresholds...')
adata.obs['passed_qc'] = True
for key in threshold_keys:
    adata.obs['passed_qc'] = adata.obs['passed_qc'] & adata.obs[key].between(*thresholds[key])

logging.info('Plot removed cells...')
plt.figure(figsize=(4, 5))
plt.grid(False)
sns.countplot(
    x='passed_qc',
    order=[True, False],
    data=adata.obs,
    hue='passed_qc',
    hue_order=[True, False],
    palette='Set2'
)
ax = plt.gca()
for pos in ['right', 'top']: 
    ax.spines[pos].set_visible(False)
for container in ax.containers:
    ax.bar_label(container)
plt.xlabel('Cell QC Status')
plt.ylabel('Count')
plt.title(f'Counts of cells QC\'d\n{dataset}')
plt.tight_layout()
plt.savefig(output_plots / 'cells_passed_all.png', bbox_inches='tight')
plt.close()


# Plot composition of removed cells
for group in groups:
    n_groups = adata.obs[group].nunique()
    if n_groups > 100:
        logging.info(f'Group {group} has too many unique values, skipping...')
        continue
    
    grouped_frac = get_fraction_removed(adata.obs, group=group, key='passed_qc')
    order = grouped_frac.sort_values('fraction_removed', ascending=False).index
    # order = adata.obs[group].value_counts().index

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(11 * (1 + n_groups/100), 6))
    sns.countplot(
        data=adata.obs,
        y=group,
        hue='passed_qc',
        hue_order=[True, False],
        order=order,
        palette='Set2',
        dodge=False,
        ax=ax1,
    )
    if n_groups < 50:
        for container in ax1.containers:
            ax1.bar_label(container)

    sns.barplot(
        data=grouped_frac,
        x='fraction_removed',
        y=group,
        order=order,
        ax=ax2,
    )
    if n_groups < 50:
        for container in ax2.containers:
            ax2.bar_label(container)
    
    for pos in ['right', 'top']: 
        ax1.spines[pos].set_visible(False)
        ax2.spines[pos].set_visible(False)
    f.suptitle(f'Cells that passed QC\n{dataset}')
    f.tight_layout()
    f.savefig(output_plots / f'by={group}.png', bbox_inches='tight')
    plt.close()

# Violin plots per QC metric
n_cols = len(threshold_keys)
fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 6))
plt.grid(False)

for i, qc_metric in enumerate(threshold_keys):
    sns.violinplot(
        data=adata.obs,
        x='passed_qc',
        order=[True, False],
        y=qc_metric,
        hue='passed_qc',
        hue_order=[True, False],
        palette='Set2',
        inner='quartile',
        legend=False,
        ax=axes[i]
    )
    axes[i].set_xlabel('QC status')
    axes[i].set_ylabel(qc_metric)
    axes[i].set_title(f'{qc_metric} distribution')
    for pos in ['right', 'top']: 
        axes[i].spines[pos].set_visible(False) 

plt.suptitle(f'Cells that passed QC\n{dataset}', fontsize=12)
plt.tight_layout()
plt.savefig(output_plots / 'per_metric_violin.png', bbox_inches='tight')
