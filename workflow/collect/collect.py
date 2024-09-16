import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from collect_utils import get_same_columns, merge_df, ALL_SLOTS


dataset = snakemake.wildcards.dataset
files = snakemake.input
output_file = snakemake.output.zarr
sep = snakemake.params.get('sep', '--')
same_slots = snakemake.params.get('same_slots', [])
merge_slots = snakemake.params.get('merge_slots', [])
kwargs = {k: k for k in same_slots+merge_slots}


if len(files) == 1:
    logging.info('Single file, write unmodified...')
    adata = read_anndata(files[0], **kwargs, dask=True, backed=True)
    print(adata, flush=True)
    if files[0].endswith('.zarr'):
        write_zarr_linked(adata, in_dir=files[0], out_dir=output_file)
    else:
        adata.write_zarr(output_file)
    exit(0)

logging.info('Read AnnData objects...')
adatas = {
    file_id: read_anndata(file, **kwargs, dask=True, backed=True)
    for file_id, file in files.items()
}

if 'obs' in merge_slots:
    same_obs_columns = get_same_columns(adatas)
else:
    same_obs_columns = []
print('same slots:', same_obs_columns)

# TODO: link slots with new slot names for merge_slots

slots = dict()
file_to_link = None

for file_id, _ad in adatas.items():

    file_name = files[file_id]
    if file_name.endswith('.zarr'):
        file_to_link = file_name

    for slot_name in merge_slots:
        slot = _ad.__dict__.get(f'_{slot_name}')
        if isinstance(slot, pd.DataFrame):
            slots[slot_name] = merge_df(
                df_current=slot,
                file_id=file_id,
                df_previous=slots.get(slot_name),
                same_columns=same_obs_columns,
                sep=sep,
            )
        elif hasattr(slot, 'items'):
            new_slot = {
                f'{key}{sep}{file_id}': value
                for key, value in slot.items()
            }
            slots[slot_name] = slots.get(slot_name, {}) | new_slot
        elif slot_name in ['X', 'raw']:
            new_slot = {
                f'{slot_name}{sep}{file_id}': _ad.__dict__.get(f'_{slot_name}')
            }
            slots['layers'] = slots.get('layers', {}) | new_slot
        else:
            raise NotImplementedError(f'Slot "{slot}" not supported')

# deal with same slots
files_to_keep = []
if file_to_link:
    files_to_keep = [
        slot for slot in ALL_SLOTS
        if slot not in same_slots
    ]
else:
    for slot_name in same_slots:
        slots[slot_name] = _ad.__dict__.get(f'_{slot_name}')

logging.info('Create AnnData object...')
adata = ad.AnnData(**slots)
print(adata, flush=True)

logging.info('Compute matrix...')
adata = dask_compute(adata)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=file_to_link,
    out_dir=output_file,
    files_to_keep=files_to_keep,
)