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
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())

keep_indices = (x for x in adata.obs_names)
remove_indices = []

if 'remove_by_column' in params:
    logging.info('remove by columns...')
    ex_filters = params['remove_by_column']
    logging.info(pformat(ex_filters))
    for column, values in ex_filters.items():
        logging.info(f'remove {values} from column="{column}"...')
        values = [str(v) for v in values]
        remove_indices.extend(
            adata.obs_names[adata.obs[column].astype(str).isin(values)]
        )

logging.info('Add filtered column...')
keep_indices = (i for i in keep_indices if i not in remove_indices)
adata.obs['filtered'] = adata.obs_names.isin(keep_indices)
value_counts = adata.obs['filtered'].value_counts()
logging.info(value_counts)

if subset and False in value_counts.index:
    kwargs |= {
        'X': 'X',
        'layers': 'layers',
        'raw': 'raw',
        'obsm': 'obsm',
        'obsp': 'obsp',
    }
    # filter out slots that aren't present in the input
    kwargs = {k: v for k, v in kwargs.items() if k in [f.name for f in Path(input_file).iterdir()]}
    
    logging.info('Read all slots for subsetting...')
    obs = adata.obs # save updated obs
    adata = read_anndata(
        input_file,
        backed=backed,
        dask=dask,
        **{k: v for k, v in kwargs.items() if k != 'obs'},
    )
    adata.obs = obs # updated obs
    
    logging.info('Subset data by filters...')
    adata = adata[adata.obs['filtered']].copy()
    logging.info(adata.__str__())

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=kwargs.keys(),
)