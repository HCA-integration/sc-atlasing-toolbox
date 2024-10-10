import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc
import scipy
import seaborn as sns
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
histplot_path = snakemake.output.histplot_path
method_name = snakemake.wildcards.method

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_zarr,
    obs='obs',
    obsp='obsp',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

dist_max = adata.obsp["distances"]
if scipy.sparse.issparse(dist_max):
    dist_max = dist_max.toarray()
nnz = (dist_max > 0).sum()
nnz_percent = round(100 * nnz / dist_max.size, 2)
sns.histplot(dist_max.flatten())
plt.title(f'{method_name}. Non zero elements: {nnz} ({nnz_percent}%)')

plt.savefig(histplot_path, bbox_inches='tight', dpi=300)
