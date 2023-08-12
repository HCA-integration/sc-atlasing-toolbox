import warnings
warnings.filterwarnings("ignore", message="Warning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params

if params is None:
    params = {}
else:
    params = {k: v for k,v in params.items()}

print('params:', params)

assert 'basis' in params, '"basis" is positional and must be defined'
basis = params['basis']

adata = read_anndata(input_file)


# convert to dense matrix
try:
    adata.obsm[basis] = np.array(adata.obsm[basis].todense())
except:
    pass

# parse colors
if 'color' in params and params['color'] is not None:
    colors = params['color'] if isinstance(params['color'], list) else [params['color']]
    # filter colors with too many categories
    params['color'] = [color for color in colors if adata.obs[color].nunique() <= 128]
    if len(params['color']) == 0:
        params['color'] = None

# plot embedding
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.embedding(adata, show=False, **params)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
