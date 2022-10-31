import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc

sns.set_theme(style='white')


def plot_qc_joint(
        adata,
        x,
        y,
        log=1,
        hue=None,
        marginal_hue=None,
        marginal_legend=False,
        palette=None,
        x_threshold=(0, np.inf),
        y_threshold=(0, np.inf),
        title=''
):
    """
    Plot scatter plot with marginal histograms from obs columns in anndata object.

    :param adata: anndata object
    :param x: obs column for x axis
    :param y: obs column for y axis
    :param log: log base for transforming values. Default 1, no transformation
    :param hue: obs column with annotations for color coding scatter plot points
    :param marginal_hue: obs column with annotations for color coding marginal plot distributions
    :param palette: a matplotlib colormap for scatterplot points
    :param x_threshold: tuple of upper and lower filter thresholds for x axis
    :param y_threshold: tuple of upper and lower filter thresholds for y axis
    :param title: Title text for plot
    """

    adata = adata.copy()

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log > 1:
        x_log = f'log1p_{x}'
        y_log = f'log1p_{y}'
        adata.obs[x_log] = log1p_base(adata.obs[x], log)
        adata.obs[y_log] = log1p_base(adata.obs[y], log)
        x_threshold = log1p_base(x_threshold, log)
        y_threshold = log1p_base(y_threshold, log)
        x = x_log
        y = y_log

    print(x)
    g = sns.JointGrid(
        data=adata.obs,
        x=x,
        y=y,
        xlim=(0, adata.obs[x].max()),
        ylim=(0, adata.obs[y].max()),
    )
    # main plot
    g.plot_joint(
        sns.scatterplot,
        data=adata[adata.obs.sample(adata.n_obs).index].obs,
        alpha=.3,
        hue=hue,
        s=2,
        palette=palette,
    )
    # marginal hist plot
    use_marg_hue = marginal_hue is not None
    g.plot_marginals(
        sns.histplot,
        data=adata.obs,
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

    return g


sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 12, 12

input_h5ad = snakemake.input.h5ad
output_joint = snakemake.output.joint
output_violin = snakemake.output.violin
output_avg = snakemake.output.average_jitter
dataset = snakemake.wildcards.dataset

print(f'Read {input_h5ad}...')
adata = sc.read(input_h5ad)

print('Calculate QC stats...')
adata.var["mito"] = adata.var['feature_name'].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

print('Joint QC plot...')
plot_qc_joint(
    adata,
    x='total_counts',
    y='n_genes_by_counts',
    hue='pct_counts_mito',
    palette='plasma',
    marginal_hue='donor',
    title=f'Joint QC for {dataset}',
)
plt.tight_layout()
plt.savefig(output_joint)

print('Violin plots...')
fig, axes = plt.subplots(nrows=2, ncols=1)

sc.pl.violin(
    adata,
    keys='total_counts',
    groupby='donor',
    rotation=90,
    show=False,
    ax=axes[0],
)

sc.pl.violin(
    adata,
    keys='pct_counts_mito',
    groupby='donor',
    rotation=90,
    show=False,
    ax=axes[1],
)
fig.suptitle(dataset)

# g = sns.violinplot(x="donor", y="total_counts", data=adata.obs)
# g.set_xticklabels(g.get_xticklabels(), rotation=90)
# g.get_figure().savefig(output_violin)

plt.tight_layout()
plt.savefig(output_violin)

print('Average scatter plots...')
fig, axes = plt.subplots(nrows=3, ncols=1)

df = (
    adata.obs[['n_genes_by_counts', 'pct_counts_mito', 'total_counts', 'donor', 'sample']]
    .groupby(['donor', 'sample'])
    .median()
    .reset_index('donor')
    .dropna()
)

sns.stripplot(
    data=df,
    x='total_counts',
    hue='donor',
    legend=False,
    ax=axes[0],
)
axes[0].set_title('Mean genes per cell, median over sample')
axes[0].set_xlim(0, None)

sns.stripplot(
    data=df,
    x='n_genes_by_counts',
    hue='donor',
    legend=False,
    ax=axes[1],
)
axes[1].set_title('Mean genes per cell, median over sample')
axes[1].set_xlim(0, None)

sns.stripplot(
    data=df,
    x='pct_counts_mito',
    hue='donor',
    legend=False,
    ax=axes[2],
)
axes[2].set_title('Percent of mitochondrial genes per cell, median over sample')
axes[2].set_xlim(0, 100)

fig.suptitle(dataset)

plt.tight_layout()
plt.savefig(output_avg)
