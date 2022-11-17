import numpy as np
from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc

sc.set_figure_params(dpi=100, frameon=False)
input_h5ad = snakemake.input.h5ad
output_png = snakemake.output.png
bulk_by = snakemake.params['bulk_by']
color = snakemake.params['color']
colors = color if isinstance(color, list) else [color]

if 'dataset' in snakemake.wildcards.keys():
    dataset = snakemake.wildcards.dataset
elif 'organ' in snakemake.wildcards.keys():
    dataset = snakemake.wildcards.organ
else:
    raise ValueError(f'Neither of "dataset" or "organ" in wildcards: {snakemake.wildcards}')

adata = sc.read(input_h5ad)

# make sure all columns are present for bulk
adata.obs['bulk_by'] = adata.obs[bulk_by].str.cat(adata.obs[color]).astype('category')


def get_pseudobulks(adata, group_key):

    pseudobulk = {'Genes': adata.var_names.values}

    for i in adata.obs.loc[:, group_key].cat.categories:
        temp = adata.obs.loc[:, group_key] == i
        pseudobulk[i] = adata[temp].X.sum(axis=0).A1

    pseudobulk = pd.DataFrame(pseudobulk).set_index('Genes')

    return pseudobulk


pbulks_df = get_pseudobulks(adata, group_key='bulk_by')

merge_cols = ['bulk_by', bulk_by]
merge_cols.extend(colors)

obs = pd.DataFrame(pbulks_df.columns, columns=['bulk_by'])
obs = obs.merge(
    adata.obs[merge_cols].drop_duplicates(),
    on='bulk_by',
    how='left'
)

adata_bulk = anndata.AnnData(
    pbulks_df.transpose().reset_index(drop=True),
    obs=obs,
    dtype='float32'
)
sc.pp.log1p(adata_bulk)
sc.pp.highly_variable_genes(adata_bulk, n_top_genes=2000)
sc.pp.pca(adata_bulk)

# remove colors if too many entries
colors = [c for c in colors if adata_bulk.obs[c].nunique() <= 64]

n_rows = len(colors)
height = np.max([7, n_rows * 0.2 * 7])
plt.rcParams['figure.figsize'] = (10, height)
sc.pl.pca(
    adata_bulk,
    color=colors,
    title=[f'{dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
    show=False,
    ncols=1
)
plt.tight_layout()
plt.savefig(output_png, bbox_inches='tight')
