import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", message="Warning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")
from matplotlib import pyplot as plt
import scanpy as sc

from utils.io import read_anndata
from utils.misc import remove_outliers


sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)

input_file = snakemake.input[0]
output_plot = Path(snakemake.output.plot)
output_additional = Path(snakemake.output.additional_plots)
output_additional.mkdir(exist_ok=True)

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

# parse neighbors key
neighbors_key = params.get('neighbors_key', 'neighbors')
if isinstance(neighbors_key, list):
    neighbors_keys = params['neighbors_key']
    del params['neighbors_key']
    for neighbors_key in neighbors_keys:
        sc.pl.embedding(adata, f'X_umap_{neighbors_key}', show=False, neighbors_key=neighbors_key, **params)
        plt.suptitle(f'neighbors_key: {neighbors_key}')
        fig_file = output_additional / f'{neighbors_key}.png'
        plt.savefig(fig_file, bbox_inches='tight', dpi=200)
    logging.info(f'link {output_plot} to {fig_file}')
    output_plot.symlink_to(fig_file.resolve(), target_is_directory=False)
else:
    # plot UMAP
    sc.pl.umap(adata, show=False, **params)
    plt.savefig(output_plot, bbox_inches='tight', dpi=200)
