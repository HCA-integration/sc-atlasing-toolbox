from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc

input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
output_file = snakemake.output[0]
output_per_group = Path(snakemake.output.per_group)
kwargs = snakemake.params.kwargs
marker_genes = snakemake.params.marker_genes

if kwargs is None:
    kwargs = {}
else:
    kwargs = {k: v for k,v in kwargs.items()}

adata = sc.read(input_file)
group_assignment = pd.read_table(input_group_assignment, index_col=0)
print(group_assignment)

# add group assignment to adata
group_cols = ['group', 'reannotation']
adata.obs[group_cols] = group_assignment.loc[adata.obs_names, group_cols]

# match marker genes and var_names
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

if isinstance(marker_genes, dict):
    marker_genes = {k: v for k, v in marker_genes.items() if k in adata.var_names}
elif isinstance(marker_genes, list):
    marker_genes = [g for g in marker_genes if g in adata.var_names]
else:
    raise ValueError('marker_genes must be dict or list')

sc.pl.dotplot(
    adata,
    groupby='group',
    var_names=marker_genes,
    show=False,
    title=f'Markers vs groups for {snakemake.wildcards.dataset}',
    **kwargs,
)
plt.savefig(output_file, bbox_inches='tight', dpi=200)

# per group
output_per_group.mkdir(exist_ok=True)

for group in adata.obs['group'].unique():
    sc.pl.dotplot(
        adata[adata.obs['group'] == group],
        groupby='reannotation',
        var_names=marker_genes,
        show=False,
        title=f'Group: {group}',
        **kwargs,
    )
    plt.savefig(output_per_group / f'group~{group}.png', bbox_inches='tight', dpi=200)
