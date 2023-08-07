from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")

from utils.misc import remove_outliers


input_h5ad = snakemake.input.h5ad
input_coordinates = snakemake.input.coordinates
input_clusters = snakemake.input.clusters
input_coordinates = snakemake.input.coordinates
output_file = snakemake.output[0]

adata = sc.read(input_h5ad)
cluster_df = pd.read_table(input_clusters, index_col=0, dtype=str)
umap = np.load(input_coordinates)

# add UMAP coordinates to adata
adata.obsm['X_umap'] = umap
# remove outliers
adata = remove_outliers(adata, 'max')
adata = remove_outliers(adata, 'min')

# add cluster assignment to adata
cluster_cols = cluster_df.columns
adata.obs[cluster_cols] = cluster_df.loc[adata.obs_names, cluster_cols]

# plot
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.umap(
    adata,
    color=cluster_cols,
    legend_loc='on data',
    ncols=2,
)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
