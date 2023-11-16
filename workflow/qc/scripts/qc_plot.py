import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata

from utils.io import read_anndata


def plot_qc_joint(
        df,
        x,
        y,
        log=1,
        hue=None,
        marginal_hue=None,
        marginal_legend=False,
        palette=None,
        x_threshold=(0, np.inf),
        y_threshold=(0, np.inf),
        title='',
        return_df=False
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log: log base for transforming values. Default 1, no transformation
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param palette: a matplotlib colormap for scatterplot points
    :param x_threshold: tuple of upper and lower filter thresholds for x axis
    :param y_threshold: tuple of upper and lower filter thresholds for y axis
    :param title: Title text for plot
    :return:
        seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log > 1:
        x_log = f'log{log} {x}'
        y_log = f'log{log} {y}'
        df[x_log] = log1p_base(df[x], log)
        df[y_log] = log1p_base(df[y], log)
        x_threshold = log1p_base(x_threshold, log)
        y_threshold = log1p_base(y_threshold, log)
        x = x_log
        y = y_log

    g = sns.JointGrid(
        data=df,
        x=x,
        y=y,
        xlim=(0, df[x].max()),
        ylim=(0, df[y].max()),
    )
    # main plot
    g.plot_joint(
        sns.scatterplot,
        data=df,
        alpha=.3,
        hue=hue,
        s=2,
        palette=palette,
    )
    # marginal hist plot
    use_marg_hue = marginal_hue is not None
    g.plot_marginals(
        sns.histplot,
        data=df,
        hue=marginal_hue,
        legend=marginal_legend,
        element='step' if use_marg_hue else 'bars',
        fill=False,
        bins=100
    )

    g.fig.suptitle(title)

    # x threshold
    for t, t_def in zip(x_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axhline(y=t, color='red')
            g.ax_marg_y.axhline(y=t, color='red')

    # y threshold
    for t, t_def in zip(y_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axvline(x=t, color='red')
            g.ax_marg_x.axvline(x=t, color='red')

    if return_df:
        return g, df
    return g


input_obs = snakemake.input.obs
output_joint = snakemake.output.joint
output_joint_mito = snakemake.output.joint_mito
output_joint_log = snakemake.output.joint_log
output_violin = snakemake.output.violin
output_avg = snakemake.output.average_jitter
hue = snakemake.params.hue
sample = snakemake.params.sample
dataset = snakemake.params.dataset

print(f'Read {input_obs}...')
obs = pd.read_table(input_obs, index_col=0)

# if no cells filtered out, save empty plots
if obs.shape[0] == 0:
    plt.savefig(output_joint)
    plt.savefig(output_joint_log)
    plt.savefig(output_violin)
    plt.savefig(output_avg)
    exit()


sns.set_theme(style='white')
sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 12, 12

print('Joint QC plot...')
plot_qc_joint(
    obs,
    x='total_counts',
    y='n_genes_by_counts',
    hue='pct_counts_mito',
    palette='plasma',
    marginal_hue=hue,
    title=f'Joint QC for {dataset}',
)
plt.tight_layout()
plt.savefig(output_joint)

plot_qc_joint(
    obs,
    x='n_genes_by_counts',
    y='pct_counts_mito',
    hue='total_counts',
    palette='plasma',
    marginal_hue=hue,
    title=f'Joint QC for {dataset}',
)
plt.tight_layout()
plt.savefig(output_joint_mito)


_, obs = plot_qc_joint(
    obs,
    x='total_counts',
    y='n_genes_by_counts',
    log=10,
    hue='pct_counts_mito',
    palette='plasma',
    marginal_hue=hue,
    title=f'Joint QC for {dataset}',
    return_df=True,
)
plt.tight_layout()
plt.savefig(output_joint_log)

print('Violin plots...')
fig, axes = plt.subplots(nrows=2, ncols=1)

adata = anndata.AnnData(obs=obs)
sc.pl.violin(
    adata,
    keys='log10 total_counts',
    groupby=hue,
    rotation=90,
    show=False,
    ax=axes[0],
)

sc.pl.violin(
    adata,
    keys='pct_counts_mito',
    groupby=hue,
    rotation=90,
    show=False,
    ax=axes[1],
)
fig.suptitle(dataset)

# g = sns.violinplot(x="donor", y="total_counts", data=obs)
# g.set_xticklabels(g.get_xticklabels(), rotation=90)
# g.get_figure().savefig(output_violin)

plt.tight_layout()
plt.savefig(output_violin)

print('Average scatter plots...')
fig, axes = plt.subplots(nrows=3, ncols=1)

df = (
    obs[['log10 n_genes_by_counts', 'pct_counts_mito', 'log10 total_counts', hue, sample]]
    .groupby([hue, sample])
    .median()
    .reset_index(hue)
    .dropna()
)

sns.stripplot(
    data=df,
    x='log10 total_counts',
    hue=hue,
    legend=False,
    ax=axes[0],
)
axes[0].set_title('Mean counts per cell, median over sample')
axes[0].set_xlim(0, None)

sns.stripplot(
    data=df,
    x='log10 n_genes_by_counts',
    hue=hue,
    legend=False,
    ax=axes[1],
)
axes[1].set_title('Mean genes per cell, median over sample')
axes[1].set_xlim(0, None)

sns.stripplot(
    data=df,
    x='pct_counts_mito',
    hue=hue,
    legend=False,
    ax=axes[2],
)
axes[2].set_title('Percent of mitochondrial genes per cell, median over sample')
axes[2].set_xlim(0, 100)

fig.suptitle(dataset)

plt.tight_layout()
plt.savefig(output_avg)
