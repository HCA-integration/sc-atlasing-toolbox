from pathlib import Path
from pandas.api.types import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns
import scanpy as sc
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds, plot_qc_joint


input_zarr = snakemake.input.zarr
output_joint = snakemake.output.joint
output_density = snakemake.output.density
output_density_log = snakemake.output.log_density
# output_violin = snakemake.output.violin
# output_avg = snakemake.output.average_jitter

output_joint = Path(output_joint)
output_joint.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')

# get parameters
file_id = snakemake.wildcards.file_id
sample, dataset, hues = parse_parameters(adata, snakemake.params)

# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('Save empty plot...')
    plt.savefig(output_density)
    exit()

thresholds = get_thresholds(
    threshold_keys=['n_counts', 'n_genes', 'percent_mito'],
    autoqc_thresholds=adata.uns['scautoqc_ranges'],
    user_thresholds=snakemake.params.get('thresholds'),
    user_key=file_id,
)
logging.info(f'\n{pformat(thresholds)}')

sns.set_theme(style='white')
sc.set_figure_params(frameon=False, fontsize=10, dpi_save=300)
plt.rcParams['figure.figsize'] = 12, 12



logging.info('Joint QC plots...')

# TODO
# density plot - done
# log scale n_genes - done
# color by donor - done
# color by disease - done

# n_counts vs n_features
x = 'n_counts'
y = 'n_genes'

for hue in hues+['percent_mito']:
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
    if is_numeric_dtype(adata.obs[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        palette = None # if adata.obs[hue].nunique() > 100 else 'plasma'
        legend = adata.obs[hue].nunique() <= 20

    plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title=joint_title,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'counts_vs_genes_hue={hue}.png')

    _, adata.obs = plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        log_x=10,
        log_y=10,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title=joint_title,
        return_df=True,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'log_counts_vs_log_genes_hue={hue}.png')


# n_feature x mito frac
x = 'n_genes'
y = 'percent_mito'

for hue in hues:
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
    if is_numeric_dtype(adata.obs[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        palette = None # if adata.obs[hue].nunique() > 100 else 'plasma'
        legend = adata.obs[hue].nunique() <= 20
    
    plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        # log_x=10,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        hue=hue,
        marginal_hue=hue,
        title=joint_title,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'genes_vs_mito_hue={hue}.png')
    
    # log n_genes
    plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        log_x=10,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        hue=hue,
        marginal_hue=hue,
        title=joint_title,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'log_genes_vs_mito_hue={hue}.png')

density_data = adata.obs.sample(n=int(min(1e5, adata.n_obs)), random_state=42)
plot_qc_joint(
    density_data,
    x=x,
    y=y,
    # main_plot_function=Axes.hexbin,
    main_plot_function=sns.kdeplot,
    # log_x=10,
    x_threshold=thresholds[x],
    y_threshold=thresholds[y],
    marginal_hue=sample,
    title=joint_title,
    fill=True,
    cmap='plasma',
    alpha=.8,
)
plt.tight_layout()
plt.savefig(output_density)

plot_qc_joint(
    density_data,
    x=x,
    y=y,
    main_plot_function=sns.kdeplot,
    log_x=10,
    x_threshold=thresholds[x],
    y_threshold=thresholds[y],
    marginal_hue=sample,
    title=joint_title,
    fill=True,
    cmap='plasma',
    alpha=.8,
)
plt.tight_layout()
plt.savefig(output_density_log)


# logging.info('Violin plots...')
# fig, axes = plt.subplots(nrows=2, ncols=1)
# hue = hues[0]
# adata.obs[hue] = adata.obs[hue].astype('category')

# sc.pl.violin(
#     adata,
#     keys='log10 n_counts',
#     groupby=hue,
#     rotation=90,
#     show=False,
#     fill=False,
#     ax=axes[0],
# )

# sc.pl.violin(
#     adata,
#     keys='percent_mito',
#     groupby=hue,
#     rotation=90,
#     show=False,
#     fill=False,
#     ax=axes[1],
# )
# fig.suptitle(dataset)

# # g = sns.violinplot(x="donor", y="n_counts", data=adata.obs)
# # g.set_xticklabels(g.get_xticklabels(), rotation=90)
# # g.get_figure().savefig(output_violin)

# plt.tight_layout()
# plt.savefig(output_violin)

# logging.info('Average scatter plots...')
# fig, axes = plt.subplots(nrows=3, ncols=1)

# df = (
#     adata.obs[['log10 n_genes', 'percent_mito', 'log10 n_counts', hue, sample]]
#     .groupby([hue, sample])
#     .median()
#     .reset_index(hue)
#     .dropna()
# )

# sns.stripplot(
#     data=df,
#     x='log10 n_counts',
#     hue=hue,
#     legend=False,
#     ax=axes[0],
# )
# axes[0].set_title('Mean counts per cell, median over sample')
# axes[0].set_xlim(0, None)

# sns.stripplot(
#     data=df,
#     x='log10 n_genes',
#     hue=hue,
#     legend=False,
#     ax=axes[1],
# )
# axes[1].set_title('Mean genes per cell, median over sample')
# axes[1].set_xlim(0, None)

# sns.stripplot(
#     data=df,
#     x='percent_mito',
#     hue=hue,
#     legend=False,
#     ax=axes[2],
# )
# axes[2].set_title('Percent of mitochondrial genes per cell, median over sample')
# axes[2].set_xlim(0, 100)

# fig.suptitle(dataset)

# plt.tight_layout()
# plt.savefig(output_avg)
