"""
Data filtering
"""
from pathlib import Path
from pprint import pformat
import numpy as np
import scanpy as sc
from dask import array as da
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = dict(snakemake.params)

backed = params.get('backed', True)
dask = params.get('dask', True)
subset = params.get('subset', False)

kwargs = {'obs': 'obs'}
if subset:
    kwargs |= {
        'X': 'X',
        'layers': 'layers',
        'obsm': 'obsm',
        'obsp': 'obsp',
    }
    # filter out slots that aren't present in the input
    kwargs = {k: v for k, v in kwargs.items() if k in [f.name for f in Path(input_file).iterdir()]}

adata = read_anndata(
    input_file,
    **kwargs,
    backed=backed,
    dask=dask
)
logging.info(adata.__str__())

keep_indices = adata.obs_names.tolist()
remove_indices = []

# explicit filter
if 'remove_by_column' in params:
    logging.info('remove by columns...')
    ex_filters = params['remove_by_column']
    logging.info(pformat(ex_filters))
    for column, values in ex_filters.items():
        logging.info(f'remove {values} in column="{column}"...')
        values = [str(v) for v in values]
        remove_indices.extend(
            adata.obs_names[adata.obs[column].astype(str).isin(values)]
        )
        remove_indices = list(set(remove_indices))

keep_indices = [i for i in keep_indices if i not in remove_indices]
n_removed = adata.n_obs - len(keep_indices)
logging.info(f'keep {len(keep_indices)} cells')
logging.info(f'remove {n_removed} cells')

logging.info('Save filtered...')
adata.obs['filtered'] = adata.obs_names.isin(keep_indices)

if subset:
    logging.info('Subset data by filters...')
    adata = adata[adata.obs['filtered']].copy()

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=kwargs.keys(),
)