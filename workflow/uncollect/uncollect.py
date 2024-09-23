import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

from utils.io import read_anndata, get_file_reader, write_zarr_linked
from utils.misc import dask_compute


def check_slot_name(slot_name, new_file_id, sep):
    splits = slot_name.split(sep)
    return len(splits) == 1 or all(x in splits for x in new_file_id.split(sep))


def remove_file_id(slot_name, new_file_id, sep, verbose=False):
    # split new_file_id, in case sep is in new_file_id
    slot_name_splits = [
        x for x in slot_name.split(sep)
        if x not in new_file_id.split(sep)
    ]
    if verbose:
        print(f'{slot_name} -> {sep.join(slot_name_splits)}', flush=True)
    return sep.join(slot_name_splits)


dataset = snakemake.wildcards.dataset
input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
new_file_id = snakemake.wildcards.new_file_id
sep = snakemake.params.get('sep', '--')


logging.info(f'Read {input_file}...')
func, _ = get_file_reader(input_file)
slots_in_file = func(input_file, 'r').keys()

kwargs = dict(dask=True, backed=True)
kwargs |= {k: k for k in slots_in_file}

if input_file.endswith('.zarr'):
    large_slots = ['layers', 'obsm', 'obsp', 'uns']
    kwargs = {k: v for k, v in kwargs.items() if k not in large_slots+['X']}
    
    adata = read_anndata(input_file, **kwargs)
    default_values = [
        csr_matrix(adata.shape),
        csr_matrix((adata.n_obs, 0)),
        csr_matrix((adata.n_obs, adata.n_obs)),
        None
    ]
    for slot_name, value in zip(large_slots, default_values):
        slot_keys = func(input_file, 'r').get(slot_name, {}).keys()
        setattr(
            adata,
            f'_{slot_name}',
            {key: value for key in slot_keys}
        )
else:
    adata = read_anndata(input_file, **kwargs)

slots = dict()
slot_link_map = dict()
files_to_keep = []

for slot_name in slots_in_file:
    logging.info(f'Process slot "{slot_name}"...')
    slot = getattr(adata, slot_name)
    
    if isinstance(slot, pd.DataFrame):
        columns = [
            col for col in slot.columns 
            if check_slot_name(col, new_file_id, sep)
        ]
        slots[slot_name] = slot[columns].rename(
            columns=lambda x: remove_file_id(x, new_file_id, sep)
        )
        files_to_keep.append(slot_name)
    
    elif hasattr(slot, 'keys'):
        if input_file.endswith('.zarr'):
            for key in slot.keys():
                if not check_slot_name(key, new_file_id, sep):
                    print('Skip', key, flush=True)
                    continue
                new_key = remove_file_id(key, new_file_id, sep)
                slot_link_map |= {f'{slot_name}/{new_key}': f'{slot_name}/{key}'}
                files_to_keep.append(slot_name)
        else:
            for key in list(slot.keys()):
                new_key = remove_file_id(key, new_file_id, sep)
                value = slot.pop(key)
                if not check_slot_name(key, new_file_id, sep):
                    slot[new_key] = value
            slots[slot_name] = slot
    
    elif slot_name in ['X', 'raw']:
        if input_file.endswith('.zarr'):
            slot_link_map |= {slot_name: slot_name}
        else:
            slots[slot_name] = slot
    else:
        raise NotImplementedError(f'Slot "{slot_name}" not supported')


logging.info('Create AnnData object...')
adata = ad.AnnData(**slots)
print(adata, flush=True)

logging.info('Compute matrix...')
adata = dask_compute(adata)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
    slot_map=slot_link_map,
)