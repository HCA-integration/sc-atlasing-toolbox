import logging
logging.basicConfig(level=logging.INFO)
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")

from utils.io import read_anndata

kwargs = dict(
    frameon=False,
    size=2,
)

input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
output_png = snakemake.output[0]
group = snakemake.wildcards.group

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm')
print(adata, flush=True)

# add group assignment to adata
logging.info('Add group assignment to adata...')
group_assignment = pd.read_table(input_group_assignment, index_col=0)
group_cols = ['group', 'reannotation']
adata.obs.loc[group_assignment.index, group_cols] = group_assignment[group_cols]

# set UMAP figure parameters
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
dot_size = np.min([np.max([0.2, 120_000 / adata.n_obs]), 200])

if group == 'all':
    logging.info('Plot global UMAP...')
    sc.pl.umap(adata, color='group', size=dot_size)
    plt.savefig(output_png, bbox_inches='tight', dpi=200)
    exit()

logging.info(f'Plot UMAP for group={group}...')
reanno_to_plot = adata[adata.obs['group'] == group].obs['reannotation'].unique().tolist()
adata.obs['reannotation'] = np.where(adata.obs['group'] == group, adata.obs['reannotation'], float('nan'))
sc.pl.umap(
    adata,
    groups=reanno_to_plot,
    color='reannotation',
    title=f'Group: {group}',
    palette=sc.pl.palettes.default_20 if adata.obs['reannotation'].nunique() <= 20 else sc.pl.palettes.default_102,
    size=dot_size,
)
plt.savefig(output_png, bbox_inches='tight', dpi=200)
