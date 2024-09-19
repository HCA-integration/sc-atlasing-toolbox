import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

from utils.io import read_anndata, get_file_reader, write_zarr_linked
from utils.misc import dask_compute

ALL_SLOTS = [
    'X',
    'obs',
    'obsm',
    'obsp',
    'var',
    'varm',
    'varp',
    'uns',
    'layers',
    'raw',
]


def check_slot_name(slot_name, new_file_id, sep):
    splits = slot_name.split(sep)
    return len(splits) == 1 or new_file_id in splits


def remove_file_id(slot_name, new_file_id, sep):
    splits = slot_name.split(sep)
    return sep.join([x for x in splits if x != new_file_id])


dataset = snakemake.wildcards.dataset
input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
new_file_id = snakemake.wildcards.new_file_id
sep = snakemake.params.get('sep', '--')


def read_file(file, **kwargs):
    if file.endswith('.zarr'):
        kwargs.pop('X', None)
        kwargs.pop('layers', None)
        adata = read_anndata(file, **kwargs)
        func, _ = get_file_reader(file)
        layers_keys = func(file, 'r').get('layers', {}).keys()
        adata.layers = {key: csr_matrix(adata.shape) for key in layers_keys}
        return adata
    return read_anndata(file, **kwargs)


logging.info('Read AnnData objects...')
adata = read_file(
    input_file,
    **{x: x for x in ALL_SLOTS},
    dask=True,
    backed=True
)

slots = dict()
slot_link_map = dict()
files_to_keep = []

for slot_name in ALL_SLOTS:
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
    
    elif hasattr(slot, 'items'):
        if input_file.endswith('.zarr'):
            for key in slot.keys():
                if not check_slot_name(key, new_file_id, sep):
                    continue
                new_key = remove_file_id(key, new_file_id, sep)
                slot_link_map |= {f'{slot_name}/{new_key}': f'{slot_name}/{key}'}
                files_to_keep.append(slot_name)
        else:
            for key in list(slot.keys()):
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