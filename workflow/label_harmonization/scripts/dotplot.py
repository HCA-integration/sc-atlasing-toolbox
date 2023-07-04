import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
from typing import MutableMapping
from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc


def filter_markers(markers, _keys):
    if isinstance(markers, dict):
        markers = {
            k: [g for g in _keys if g in v]
            for k, v in markers.items()
        }
        return {k: l for k, l in markers.items() if len(l) > 0}
    elif isinstance(markers, list):
        return [g for g in markers if g in _keys]
    else:
        raise ValueError('marker must be dict or list')


def mapping_to_dict(mapping: MutableMapping) -> dict:
    return {} if mapping is None else dict(mapping.items())


input_file = snakemake.input[0]
input_group_assignment = snakemake.input.group_assignment
output_file = snakemake.output[0]
output_per_group = Path(snakemake.output.per_group)
output_per_group.mkdir(exist_ok=True)

kwargs = mapping_to_dict(snakemake.params.kwargs)
marker_genes = snakemake.params.marker_genes

adata = sc.read(input_file)
group_assignment = pd.read_table(input_group_assignment, index_col=0)

# add group assignment to adata
group_cols = ['group', 'reannotation']
adata.obs[group_cols] = group_assignment.loc[adata.obs_names, group_cols]

# match marker genes and var_names
logging.info(adata.var)
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

logging.info('Marker genes before filtering:')
logging.info(marker_genes)
marker_genes = filter_markers(marker_genes, adata.var_names)
logging.info('Marker genes after filtering:')
logging.info(marker_genes)

logging.info('Plotting...')
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
for group in adata.obs['group'].unique():
    ad = adata[adata.obs['group'] == group]
    if ad.n_obs == 0: continue

    logging.info(f'Plotting for {group}...')
    sc.pl.dotplot(
        ad,
        groupby='reannotation',
        var_names=filter_markers(marker_genes, ad.var_names),
        show=False,
        title=f'Group: {group}',
        **kwargs,
    )
    plt.savefig(output_per_group / f'group~{group}.png', bbox_inches='tight', dpi=200)
