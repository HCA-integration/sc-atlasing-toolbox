"""
Data filtering
"""
from pprint import pformat
import numpy as np
import scanpy as sc
from dask import array as da
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, link_zarr, to_memory
from utils.misc import apply_layers


def save_empty(input_zarr, output_zarr, output_removed):
    # link unfiltered
    link_zarr(
        in_dir=input_zarr,
        out_dir=output_zarr,
    )

    # save removed cells
    adata_empty = read_anndata(input_zarr, var='var', uns='uns')
    adata_empty.X = np.zeros(shape=(0, adata_empty.var.shape[0]))
    adata_empty.write_zarr(output_removed)
    exit(0)


input_zarr = snakemake.input[0]
output_zarr = snakemake.output[0]
# output_removed = snakemake.output.removed
params = dict(snakemake.params)

if len(params) == 0:
    save_empty(input_zarr, output_zarr, output_removed)

backed = params.get('backed', True)
dask = params.get('dask', True)
adata = read_anndata(input_zarr, backed=backed, dask=dask)
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
adata_removed = adata.n_obs - len(keep_indices)
logging.info(f'keep {len(keep_indices)} cells')
logging.info(f'remove {adata_removed} cells')

logging.info('Save filtered...')
rechunk = lambda x: x.rechunk({0: 1000}) if isinstance(x, da.Array) else x
adata_filtered = apply_layers(adata[adata.obs_names.isin(keep_indices)], func=rechunk)
adata_filtered.write_zarr(output_zarr)

# logging.info('Save removed...')
# if adata_removed == 0:
#     save_empty(input_zarr, output_zarr, output_removed)
# adata_removed = apply_layers(adata[~adata.obs_names.isin(keep_indices)], func=rechunk)
# adata_removed.write_zarr(output_removed)
