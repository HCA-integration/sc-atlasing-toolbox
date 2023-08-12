import sys
import warnings
warnings.filterwarnings("ignore", message="Warning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

from utils.io import read_anndata
from utils.misc import remove_outliers


input_file = snakemake.input[0]
output_plot = snakemake.output.plot
params = {k: v for k, v in snakemake.params.items()}
if 'outlier_factor' in params:
    outlier_factor = params['outlier_factor']
    del params['outlier_factor']
else:
    outlier_factor = 3

adata = read_anndata(input_file)

# remove outliers
adata = remove_outliers(adata, 'max', factor=outlier_factor)
adata = remove_outliers(adata, 'min', factor=outlier_factor)

# parse colors
if 'color' in params and params['color'] is not None:
    colors = params['color'] if isinstance(params['color'], list) else [params['color']]
    # filter colors with too many categories
    params['color'] = [color for color in colors if adata.obs[color].nunique() <= 128]
    if len(params['color']) == 0:
        params['color'] = None

# plot UMAP
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.umap(adata, show=False, **params)
plt.savefig(output_plot, bbox_inches='tight', dpi=200)
