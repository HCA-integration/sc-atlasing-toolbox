from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")

from utils.io import read_anndata
from utils.misc import remove_outliers


kwargs = dict(
    frameon=False,
    size=2,
)

input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
output_png = snakemake.output[0]
output_per_group = Path(snakemake.output.per_group)
output_per_group.mkdir(exist_ok=True)

adata = read_anndata(input_file, obs='obs', obsm='obsm', var='var')
print(adata, flush=True)

# add grouq assignment to adata
group_assignment = pd.read_table(input_group_assignment, index_col=0)
group_cols = ['group', 'reannotation']
adata.obs.loc[group_assignment.index, group_cols] = group_assignment[group_cols]

# remove outliers
adata = remove_outliers(adata, 'max')
adata = remove_outliers(adata, 'min')

# plot
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.umap(adata, color='group')
plt.savefig(output_png, bbox_inches='tight', dpi=200)

# per group
for group in adata.obs['group'].unique():
    adata.obs['reannotation_tmp'] = np.where(adata.obs['group'] == group, adata.obs['reannotation'], '')
    sc.pl.umap(
        adata,
        groups=adata[adata.obs['group'] == group].obs['reannotation'].unique(),
        color='reannotation_tmp',
        title=f'Group: {group}',
        palette=sc.pl.palettes.default_20 if adata.obs['reannotation_tmp'].nunique() <= 20 else sc.pl.palettes.default_102,
    )
    plt.savefig(output_per_group / f'group~{group}.png', bbox_inches='tight', dpi=200)
    del adata.obs['reannotation_tmp']
