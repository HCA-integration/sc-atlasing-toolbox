from matplotlib import pyplot as plt
import scanpy as sc
import anndata

sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 30, 25
input_file = snakemake.input.zarr
output_png = snakemake.output.png
dataset = snakemake.params.dataset
markers = snakemake.params.markers

adata = anndata.read_zarr(input_file)

# if no cells filtered out, save empty plots
if adata.n_obs == 0:
    plt.savefig(output_png)
    exit()


author_label = 'author_annotation'
ontology_label = 'cell_type'
print(adata)
# match marker genes and var_names
adata.var_names = adata.var['feature_name'].astype(str)
print({k: len(v) for k, v in markers.items()})
markers = {
    k: adata.var[adata.var['feature_name'].isin(v)].index.to_list()
    for k, v in markers.items()
}
print({k: len(v) for k, v in markers.items()})

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
