import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute
from subset_functions import SUBSET_MAP

input_file = snakemake.input[0]
output_file = snakemake.output[0]
strategy = snakemake.wildcards.strategy
n_cell_max = snakemake.params.get('n_cells')
n_cell_max = np.iinfo(int).max if n_cell_max is None else int(n_cell_max)
sample_key = snakemake.params.get('sample_key')
# label_key = snakemake.params.get('label_key')

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, backed=True, dask=True)
logging.info(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    logging.info(adata.__str__())
    raise AssertionError(f'sample key "{sample_key}" not in adata')

# TODO: call function
subset_func = SUBSET_MAP.get(strategy, ValueError(f'Unknown strategy: {strategy}'))
subset_mask = subset_func(adata, n_cell_max, sample_key)

# subset_key = f'subset_{strategy}'
# adata.obs[subset_key] = False
# adata.obs.loc[subset_mask, subset_key] = True

n_cells = subset_mask.sum()
n_samples = len(adata.obs.loc[subset_mask, sample_key].unique())

logging.info(f'Subset to {n_cells} cells, {n_samples} samples...')
adata = adata[subset_mask].copy()
adata.uns['subset'] = strategy
logging.info(adata.__str__())

logging.info('Compute dask array...')
dask_compute(adata)

# save
logging.info('Write...')
adata.write_zarr(output_file)
