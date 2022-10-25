from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc

sc.set_figure_params(dpi=100, frameon=False)
plt.rcParams['figure.figsize'] = 20, 12
input_h5ad = snakemake.input.h5ad
output_png = snakemake.output.png
dataset = snakemake.wildcards.dataset
markers = snakemake.params['markers']

adata = sc.read(input_h5ad)

author_label = adata.uns['meta']['cell_annotation']
ontology_label = 'cell_type'

# match marker genes and var_names
adata.var_names = adata.var['feature_name'].astype(str)
n_markers_all = len(markers)
markers = adata.var[adata.var['feature_name'].isin(markers)].index.to_list()
print(f'removed {n_markers_all - len(markers)} genes from marker list')

fig, axes = plt.subplots(nrows=2, ncols=1)
sc.pl.dotplot(
    adata,
    markers,
    groupby=author_label,
    use_raw=False,
    standard_scale='var',
    title=f'Marker genes on {dataset} for author cell types, colummn: "{author_label}"',
    show=False,
    ax=axes[0],
)
sc.pl.dotplot(
    adata,
    markers,
    groupby=ontology_label,
    use_raw=False,
    standard_scale='var',
    title=f'Marker genes on {dataset} for ontology cell types, column: "{ontology_label}"',
    show=False,
    ax=axes[1],
)

plt.tight_layout()
plt.savefig(output_png)
