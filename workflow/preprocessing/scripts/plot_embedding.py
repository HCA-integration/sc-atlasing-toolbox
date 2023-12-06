import warnings
warnings.filterwarnings("ignore", message="Warning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params

wildcards_string = ', '.join([f'{k}: {v}' for k, v in snakemake.wildcards.items()])
if params is None:
    params = {}
else:
    params = {k: v for k,v in params.items()}

print('params:', params)

assert 'basis' in params, '"basis" is positional and must be defined'
basis = params['basis']

adata = read_anndata(input_file, obs='obs', obsm='obsm')

if adata.obs.shape[0] == 0:
    plt.savefig(output_file)
    exit()

# convert to dense matrix
try:
    adata.obsm[basis] = np.array(adata.obsm[basis].todense())
except:
    pass

# parse colors
if 'color' in params and params['color'] is not None:
    colors = params['color'] if isinstance(params['color'], list) else [params['color']]
    # remove that are not in the data
    colors = [color for color in colors if color in adata.obs.columns]
    # filter colors with too few or too many categories
    colors = [color for color in colors if 1 < adata.obs[color].nunique() <= 128 or is_numeric_dtype(adata.obs[color])]
    
    # set color parameters
    if len(colors) > 0:
        for color in colors:
            if adata.obs[color].dtype.name == 'category':
                adata.obs[color] = adata.obs[color].astype('str')
        params['color'] = colors
    else:
        del params['color']

# plot embedding
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.embedding(adata, show=False, **params)
plt.suptitle(f'{wildcards_string}, n={adata.n_obs}')
plt.savefig(output_file, bbox_inches='tight', dpi=200)
