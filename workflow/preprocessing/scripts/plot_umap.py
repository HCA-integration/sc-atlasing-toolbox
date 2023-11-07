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

wildcards_string = ', '.join([f'{k}: {v}' for k, v in snakemake.wildcards.items()])
params = dict(snakemake.params.items())
if 'outlier_factor' in params:
    outlier_factor = params['outlier_factor']
    del params['outlier_factor']
else:
    outlier_factor = 10

adata = read_anndata(input_file, obs='obs', obsm='obsm')

# parse colors
if 'color' in params and params['color'] is not None:
    colors = params['color'] if isinstance(params['color'], list) else [params['color']]
    # remove that are not in the data
    colors = [color for color in colors if color in adata.obs.columns]
    # filter colors with too few or too many categories
    params['color'] = [color for color in colors if 1 < adata.obs[color].nunique() <= 128]
    if len(params['color']) == 0:
        params['color'] = None
    else:
        for color in params['color']:
            if adata.obs[color].dtype.name == 'category':
                adata.obs[color] = adata.obs[color].astype('str')

# parse neighbors key
neighbors_key = params.get('neighbors_key', 'neighbors')
if isinstance(neighbors_key, list):
    neighbors_keys = params['neighbors_key']
    del params['neighbors_key']
    for neighbors_key in neighbors_keys:
        basis = f'X_umap_{neighbors_key}'
        # remove outliers
        adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
        adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)
        sc.pl.embedding(
            adata[adata.obs.sample(adata.n_obs).index],
            basis,
            show=False,
            neighbors_key=neighbors_key,
            **params
        )
        plt.suptitle(f'{wildcards_string}, neighbors_key: {neighbors_key}, n={adata.n_obs}')
        fig_file = output_additional / f'{neighbors_key}.png'
        plt.savefig(fig_file, bbox_inches='tight', dpi=200)
    logging.info(f'link {output_plot} to {fig_file}')
    output_plot.symlink_to(fig_file.resolve(), target_is_directory=False)
else:
    # plot UMAP
    basis = 'X_umap'
    # remove outliers
    adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
    adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)
    sc.pl.umap(
        adata[adata.obs.sample(adata.n_obs).index],
        show=False,
        **params
    )
    plt.suptitle(f'{wildcards_string}, n={adata.n_obs}')
    plt.savefig(output_plot, bbox_inches='tight', dpi=200)
