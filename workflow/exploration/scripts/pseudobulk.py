import numpy as np
from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc

import warnings
warnings.simplefilter("ignore", UserWarning)

def get_pseudobulks(adata, group_key, agg='mean'):
    pseudobulk = {'Genes': adata.var_names.values}

    for i in adata.obs.loc[:, group_key].cat.categories:
        temp = adata.obs.loc[:, group_key] == i
        if agg == 'sum':
            pseudobulk[i] = adata[temp].X.sum(axis=0).A1
        elif agg == 'mean':
            pseudobulk[i] = adata[temp].X.mean(axis=0).A1
        else:
            raise ValueError(f'invalid aggregation method "{agg}"')

    pseudobulk = pd.DataFrame(pseudobulk).set_index('Genes')

    return pseudobulk


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_pca_1_2 = snakemake.output.pca_1_2
output_pca_2_3 = snakemake.output.pca_2_3
output_pca_scree = snakemake.output.pca_scree
bulk_by = snakemake.params['bulk_by']
color = snakemake.params['color']
colors = color if isinstance(color, list) else [color]
dataset = snakemake.params.dataset

adata = anndata.read_zarr(input_zarr)

if adata.n_obs == 0:
    plt.savefig(output_pca_1_2)
    plt.savefig(output_pca_2_3)
    plt.savefig(output_pca_scree)
    exit(0)

# normalize counts
sc.pp.normalize_total(adata)

# make sure all columns are present for bulk
adata.obs['bulk_by'] = adata.obs[bulk_by].str.cat(adata.obs[color].astype(str)).astype('category')

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

# process pseudobulk adata
sc.pp.log1p(adata_bulk)
sc.pp.highly_variable_genes(adata_bulk, n_top_genes=2000)
sc.pp.pca(adata_bulk)

# scree plo
plt.plot(adata_bulk.uns['pca']['variance_ratio'], marker='.')
plt.title('Scree plot: PC variance contribution')
plt.savefig(output_pca_scree)

# remove empty columns
colors = [c for c in colors if not adata_bulk.obs[c].isna().all()]
# remove colors if too many entries
colors = [c for c in colors if adata_bulk.obs[c].nunique() <= 64]

n_rows = len(colors)
height = np.max([7, n_rows * 0.2 * 7])
plt.rcParams['figure.figsize'] = (10, height)

sc.pl.pca(
    adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
    components='1,2',
    color=colors,
    title=[f'PCA 1-2 {dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
    show=False,
    ncols=1
)
plt.tight_layout()
plt.savefig(output_pca_1_2, bbox_inches='tight')

sc.pl.pca(
    adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
    components='2,3',
    color=colors,
    title=[f'PCA 2-3{dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
    show=False,
    ncols=1
)
plt.tight_layout()
plt.savefig(output_pca_2_3, bbox_inches='tight')
