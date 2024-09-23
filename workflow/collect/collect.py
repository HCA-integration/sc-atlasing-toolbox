import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
from pprint import pformat
from scipy.sparse import csr_matrix

from utils.io import read_anndata, get_file_reader, write_zarr_linked
from utils.misc import dask_compute
from collect_utils import get_same_columns, merge_df


dataset = snakemake.wildcards.dataset
files = snakemake.input
output_file = snakemake.output.zarr
sep = snakemake.params.get('sep', '--')
obs_index_col_map = snakemake.params.get('obs_index_col', {})
same_slots = snakemake.params.get('same_slots', [])
merge_slots = snakemake.params.get('merge_slots', [])
kwargs = {k: k for k in same_slots+merge_slots}

# parse obs index
if isinstance(obs_index_col_map, str):
    obs_index_col_map = {file_id: obs_index_col_map for file_id in files}
assert isinstance(obs_index_col_map, dict), 'obs_index_col_map must be a dict'

if len(files) == 1:
    logging.info('Single file, write unmodified...')
    adata = read_anndata(files[0], **kwargs, dask=True, backed=True)
    print(adata, flush=True)
    if files[0].endswith('.zarr'):
        write_zarr_linked(adata, in_dir=files[0], out_dir=output_file)
    else:
        adata.write_zarr(output_file)
    exit(0)


def read_file(file, **kwargs):
    logging.info(f'Read {file}...')
    if file.endswith('.zarr'):
        large_slots = ['layers', 'obsm', 'obsp', 'uns']
        kwargs = {k: v for k, v in kwargs.items() if k not in large_slots+['X']}
        
        adata = read_anndata(file, **kwargs)
        func, _ = get_file_reader(file)
        
        default_values = [
            csr_matrix(adata.shape),
            csr_matrix((adata.n_obs, 0)),
            csr_matrix((adata.n_obs, adata.n_obs)),
            None
        ]
        for slot_name, value in zip(large_slots, default_values):
            slot_keys = func(file, 'r').get(slot_name, {}).keys()
            setattr(
                adata,
                f'_{slot_name}',
                {key: value for key in slot_keys}
            )
        
        return adata
    return read_anndata(file, **kwargs)


adatas = {
    file_id: read_file(file, **kwargs, dask=True, backed=True)
    for file_id, file in files.items()
}

if 'obs' in merge_slots:
    for file_id, _ad in adatas.items():
        obs_index_column = obs_index_col_map.get(file_id)
        if obs_index_column is not None:
            logging.info(f'Set obs index "{obs_index_column}" for file_id={file_id}...')
            assert obs_index_column in _ad.obs.columns, \
                f'Index column "{obs_index_column}" not found for {file_id}\n{_ad.obs}'
            adatas[file_id].obs = _ad.obs.set_index(obs_index_column)
    logging.info('Determine which columns are the same...')
    same_obs_columns = get_same_columns(adatas, n_threads=snakemake.threads)
    logging.info(f'Same columns:\n{pformat(same_obs_columns)}')
else:
    same_obs_columns = []

# TODO: link slots with new slot names for merge_slots

slots = dict()
slot_link_map = dict()
in_dir_map = dict()
file_to_link = None

for file_id, _ad in adatas.items():

    file_name = files[file_id]
    # intialize file_to_link
    if not file_to_link and file_name.endswith('.zarr'):
        file_to_link = file_name

    for slot_name in merge_slots:
        slot = _ad.__dict__.get(f'_{slot_name}')
        update_slot_link_map = dict()

        if isinstance(slot, pd.DataFrame):
            slots[slot_name] = merge_df(
                df_current=slot,
                file_id=file_id,
                df_previous=slots.get(slot_name),
                same_columns=same_obs_columns,
                sep=sep,
            )
        elif slot_name == 'X':
            if file_name.endswith('.zarr'):
                update_slot_link_map[f'layers/{slot_name}{sep}{file_id}'] = f'{slot_name}'
            else:
                new_slot = {f'{slot_name}{sep}{file_id}': slot}
                slots['layers'] = slots.get('layers', {}) | new_slot
        elif hasattr(slot, 'items'):
            if file_name.endswith('.zarr'):
                update_slot_link_map = {
                    f'{slot_name}/{key}{sep}{file_id}': f'{slot_name}/{key}'
                    for key in slot.keys()
                }
            else:
                new_slot = {
                    f'{key}{sep}{file_id}': value
                    for key, value in slot.items()
                }
                slots[slot_name] = slots.get(slot_name, {}) | new_slot
        else:
            raise NotImplementedError(f'Slot "{slot}" not supported')
        
        slot_link_map |= update_slot_link_map
        in_dir_map |= {
            f'{slot_name}/{key}': file_name for key in update_slot_link_map.keys()
        }

# deal with same slots
files_to_keep = []
if file_to_link:
    files_to_keep = [slot for slot in slots]
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
    slot_map=slot_link_map,
    in_dir_map=in_dir_map,
)