"""
Data filtering
"""
from pprint import pformat
import numpy as np
import scanpy as sc
from dask import array as da
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.misc import apply_layers


input_file = snakemake.input[0]
output_file = snakemake.output[0]
# output_removed = snakemake.output.removed
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
            adata[adata.obs[column].astype(str).isin(values)].obs_names
        )
        remove_indices = list(set(remove_indices))

# # remove sample with different number of cells per sample
# if 'cells_per_sample' in params.keys():
#     min_cells = int(params['cells_per_sample'].get('min', 0))
#     max_cells = int(params['cells_per_sample'].get('max', np.inf))
#     logging.info(f'remove by min={min_cells}, max={max_cells} cells per sample...')
#     n_cells_per_sample = adata.obs['sample'].value_counts()
#     samples_to_keep = n_cells_per_sample[n_cells_per_sample.between(min_cells, max_cells)]
#     keep_indices.extend(
#         adata[adata.obs['sample'].isin(samples_to_keep)].obs_names
#     )

keep_indices = [i for i in keep_indices if i not in remove_indices]
n_removed = adata.n_obs - len(keep_indices)
logging.info(f'keep {len(keep_indices)} cells')
logging.info(f'remove {n_removed} cells')

logging.info('Save filtered...')
rechunk = lambda x: x.rechunk({0: 1000}) if isinstance(x, da.Array) else x
adata.obs['filtered'] = adata.obs_names.isin(keep_indices)

if subset:
    logging.info('Subset data by filters...')
    adata = adata[adata.obs['filtered']].copy()
    adata = apply_layers(adata, func=rechunk)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=kwargs.keys(),
)

# logging.info('Save removed...')
# if n_removed =(input= 0:
#     save_empty_zarr, output_file, output_removed)
# adata_removed = apply_layers(adata[~adata.obs_names.isin(keep_indices)], func=rechunk)
# adata_removed.write_zarr(output_removed)
