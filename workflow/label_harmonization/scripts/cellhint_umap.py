import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")

from utils.io import read_anndata


input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
output_dir = snakemake.output[0]
group = snakemake.wildcards.group
author_label_key = snakemake.params.author_label_key
dataset_key = snakemake.params.dataset_key

output_dir = Path(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm')
print(adata, flush=True)

# add group assignment to adata
logging.info('Add group assignment to adata...')
group_assignment = pd.read_table(input_group_assignment, index_col=0)
group_cols = ['group', 'reannotation']
adata.obs.loc[group_assignment.index, group_cols] = group_assignment[group_cols]

# shuffle cells
adata = adata[np.random.permutation(adata.n_obs)].copy()

# set UMAP figure parameters
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
dot_size = min(max(0.2, 500_000 / adata.n_obs), 200)

kwargs = dict(
    frameon=False,
    size=dot_size,
    sort_order=True,
)

png_kwargs = dict(
    dpi=300,
    bbox_inches='tight'
)

if group == 'all':
    logging.info('Plot global UMAP...')

    sc.pl.umap(adata, color='group', **kwargs)
    plt.savefig(output_dir / 'group.png', **png_kwargs)
    
    sc.pl.umap(adata, color='reannotation', **kwargs)
    plt.savefig(output_dir / 'reannotation.png', **png_kwargs)

    sc.pl.umap(adata, color=dataset_key, **kwargs | {'sort_order': False})
    plt.savefig(output_dir / 'dataset.png', **png_kwargs)

    for dataset in adata.obs[dataset_key].unique():
        groups = adata.obs[adata.obs[dataset_key] == dataset][author_label_key].unique()
        sc.pl.umap(adata, color=author_label_key, groups=groups, **kwargs)
        plt.savefig(output_dir / f'dataset={dataset}.png', **png_kwargs)

    exit()

logging.info(f'Plot UMAP for group={group}...')
# reanno_to_plot = adata[adata.obs['group'] == group].obs['reannotation'].unique().tolist()
adata.obs[author_label_key] = adata.obs[author_label_key].astype(str)
adata.obs.loc[adata.obs['group'] != group, ['reannotation', author_label_key]] = float('nan')
palette = sc.pl.palettes.default_20 if adata.obs['reannotation'].nunique() <= 20 else sc.pl.palettes.default_102

sc.pl.umap(
    adata,
    # groups=reanno_to_plot,
    color='reannotation',
    title=f'Group: {group}',
    palette=palette,
    **kwargs
)
plt.savefig(output_dir / 'reannotation.png', **png_kwargs)

sc.pl.umap(
    adata,
    # groups=reanno_to_plot,
    color=author_label_key,
    title=f'Group: {group}',
    palette=palette,
    **kwargs
)
plt.savefig(output_dir / f'{author_label_key}.png', **png_kwargs)
