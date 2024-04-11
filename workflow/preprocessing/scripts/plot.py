import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import traceback
import warnings
warnings.filterwarnings("ignore")

from matplotlib import pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype, is_string_dtype, is_categorical_dtype
from pprint import pformat

from utils.io import read_anndata
from utils.misc import ensure_dense, remove_outliers


sc.set_figure_params(
    frameon=False,
    vector_friendly=True,
    fontsize=9,
    figsize=(6,6),
    dpi=300,
    dpi_save=300
)

input_file = snakemake.input[0]
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(exist_ok=True)

params = dict(snakemake.params.items())
basis = params['basis']
wildcards_string = '\n'.join([f'{k}: {v}' for k, v in snakemake.wildcards.items()])
logging.info(f'Wildcard string: {wildcards_string}...')

logging.info(f'Read file {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm')
ensure_dense(adata, basis)

if adata.n_obs == 0:
    logging.info('No cells, skip...')
    exit()

# parse colors
colors = params.get('color', [None])
if 'color' in params:
    logging.info(f'Configured colors:\n{pformat(colors)}')
    colors = colors if isinstance(colors, list) else [colors]
    # remove that are not in the data
    colors = [color for color in colors if color in adata.obs.columns]
    # filter colors with too few or too many categories
    logging.info(f'Colors after filtering:\n{pformat(colors)}')
    
    # else:
    if len(colors) == 0:
        logging.info('No valid colors, skip...')
        colors = [None]
    else:
        for color in colors:
            column = adata.obs[color]
            if is_categorical_dtype(column) or is_string_dtype(column):
                column = column.replace(['NaN', 'None', '', 'nan'], float('nan'))
                column = pd.Categorical(column)
                adata.obs[color] = column.codes if len(column.categories) > 102 else column
    del params['color']

# remove outliers
outlier_factor = params.get('outlier_factor', 0)
params.pop('outlier_factor', None)

adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)

# set minimum point size
default_size = 150000 / adata.n_obs
size = params.get('size', default_size)
if size is None:
    size = default_size
params['size'] = np.min([np.max([size, 3, default_size]), 200])
print(f'Size: {params["size"]}', flush=True)
print(default_size, flush=True)
params['size'] = None


for color in colors:
    logging.info(f'Plot color "{color}"...')
    palette = sc.pl.palettes.godsnot_102 if is_categorical_dtype(adata.obs[color]) and adata.obs[color].nunique() > 20 else None
    try:
        fig = sc.pl.embedding(
            adata[adata.obs.sample(adata.n_obs).index],
            color=color,
            show=False,
            return_fig=True,
            palette=palette,
            **params
        )
        fig.suptitle(f'{wildcards_string}\nn={adata.n_obs}')
        legend = fig.get_axes()[0].get_legend()
        if legend:
            legend_bbox = legend.get_window_extent()
            fig_width, fig_height = fig.get_size_inches()
            fig_width = fig_width + (legend_bbox.width / fig.dpi)
            fig.set_size_inches((fig_width, fig_height))
        fig.tight_layout()
    except Exception as e:
        logging.error(f'Failed to plot {color}: {e}')
        traceback.print_exc()
        plt.plot([])
    
    out_path = output_dir / f'{color}.png'
    plt.savefig(out_path, bbox_inches='tight')
