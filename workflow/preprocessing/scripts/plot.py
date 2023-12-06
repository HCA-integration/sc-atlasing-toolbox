import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", message="Warning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")
from matplotlib import pyplot as plt
import scanpy as sc
from pandas.api.types import is_numeric_dtype

from utils.io import read_anndata
from utils.misc import ensure_dense, remove_outliers


sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)

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

if adata.obs.shape[0] == 0:
    logging.info('No cells, skip...')
    exit()

# parse colors
colors = params.get('color', [None])
if 'color' in params:
    colors = colors if isinstance(colors, list) else [colors]
    # remove that are not in the data
    colors = [color for color in colors if color in adata.obs.columns]
    # filter colors with too few or too many categories
    colors = [color for color in colors if 1 < adata.obs[color].nunique() <= 128 or is_numeric_dtype(adata.obs[color])]
    
    # set color parameters
    if len(colors) > 0:
        for color in colors:
            if adata.obs[color].dtype.name == 'category':
                adata.obs[color] = adata.obs[color].astype('str')
    else:
        colors = [None]
    del params['color']

# remove outliers
if 'outlier_factor' in params:
    outlier_factor = params['outlier_factor']
    del params['outlier_factor']
else:
    outlier_factor = 0
adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)

for color in colors:
    logging.info(f'Plot color {color}...')
    sc.pl.embedding(
        adata[adata.obs.sample(adata.n_obs).index],
        color=color,
        show=False,
        **params
    )
    plt.suptitle(f'{wildcards_string}\nn={adata.n_obs}')
    plt.tight_layout()
    
    out_path = output_dir / f'{color}.png'
    plt.savefig(out_path, bbox_inches='tight', dpi=200)
