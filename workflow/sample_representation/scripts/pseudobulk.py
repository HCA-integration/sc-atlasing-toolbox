import numpy as np
from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


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
output_h5ad = snakemake.output.h5ad
bulk_by = snakemake.params.get('bulk_by')
dataset = snakemake.wildcards.file_id

logging.info(f'Read "{input_zarr}"...')
adata = read_anndata(input_zarr, X='X', obs='obs', var='var')

# normalize counts
# sc.pp.normalize_total(adata)

# make sure all columns are present for bulk
logging.info(f"Number of pseudobulks: {adata.obs[bulk_by].nunique()}")

pbulks_df = get_pseudobulks(adata, group_key=bulk_by)

obs = pd.DataFrame(pbulks_df.columns, columns=[bulk_by]).merge(
    adata.obs.drop_duplicates(subset=[bulk_by]),
    on=bulk_by,
    how='right'
)
logging.info('Reset categories...')
for col in obs.columns:
    if obs[col].dtype.name == 'category':
        obs[col] = obs[col].astype(str).astype('category')

adata_bulk = anndata.AnnData(
    pbulks_df.transpose().reset_index(drop=True),
    obs=obs,
    dtype='float32'
)

logging.info(f'Write "{output_h5ad}"...')
adata_bulk.write(output_h5ad)

# # process pseudobulk adata
# sc.pp.log1p(adata_bulk)
# sc.pp.highly_variable_genes(adata_bulk, n_top_genes=2000)
# sc.pp.pca(adata_bulk)

# # scree plot TODO: move to preprocessing plots
# plt.plot(adata_bulk.uns['pca']['variance_ratio'], marker='.')
# plt.title('Scree plot: PC variance contribution')
# plt.savefig(output_pca_scree)

# # Check if it's enough for at least the first two componenets
# if adata.obs[bulk_by].nunique() < 3:
#     plt.savefig(output_pca_1_2)
#     plt.savefig(output_pca_2_3)
#     exit(0)

# # remove empty columns
# colors = [c for c in colors if not adata_bulk.obs[c].isna().all()]
# # remove colors if zero or too many entries
# colors = [c for c in colors if 0 < adata_bulk.obs[c].nunique() <= 64]
# # plot without colors, if no colors are available
# if len(colors) == 0:
#     colors = None

# n_rows = len(colors)
# height = np.max([7, n_rows * 0.2 * 7])
# plt.rcParams['figure.figsize'] = (10, height)

# sc.pl.pca(
#     adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
#     components='1,2',
#     color=colors,
#     title=[f'PCA 1-2 {dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
#     show=False,
#     ncols=1
# )
# plt.tight_layout()
# plt.savefig(output_pca_1_2, bbox_inches='tight')

# # missing third component
# if len(adata_bulk.uns['pca']['variance_ratio']) < 3:
#     plt.savefig(output_pca_2_3)
#     exit(0)

# sc.pl.pca(
#     adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
#     components='2,3',
#     color=colors,
#     title=[f'PCA 2-3{dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
#     show=False,
#     ncols=1
# )
# plt.tight_layout()
# plt.savefig(output_pca_2_3, bbox_inches='tight')
