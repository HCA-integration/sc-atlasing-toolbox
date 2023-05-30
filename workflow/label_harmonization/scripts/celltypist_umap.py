from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

kwargs = dict(
    frameon=False,
    size=2,
)

input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
input_coordinates = snakemake.input.coordinates
output_file = snakemake.output[0]
output_per_group = Path(snakemake.output.per_group)
output_per_group.mkdir(exist_ok=True)

adata = sc.read(input_file)
group_assignment = pd.read_table(input_group_assignment, index_col=0)
umap = np.load(input_coordinates)

# add UMAP coordinates to adata
adata.obsm['X_umap'] = umap

# add grouq assignment to adata
group_cols = ['group', 'reannotation']
adata.obs[group_cols] = group_assignment.loc[adata.obs_names, group_cols]

sc.pl.umap(adata, color='group')
plt.savefig(output_file, bbox_inches='tight', dpi=200)

# per group
for group in adata.obs['group'].unique():
    sc.pl.umap(
        adata[adata.obs['group'] == group],
        color='reannotation',
    )
    plt.savefig(output_per_group / f'group~{group}.png', bbox_inches='tight', dpi=200)
