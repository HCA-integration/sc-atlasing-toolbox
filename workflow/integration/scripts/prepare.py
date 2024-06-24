from anndata import AnnData
from dask import array as da
import numpy as np
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute
from utils.processing import assert_neighbors, sc


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
batch_keys = params.batches
label_keys = params.labels
var_mask = snakemake.wildcards.var_mask
save_subset = params.get('save_subset', False)
recompute_pca = params.get('recompute_pca', True)


def read_and_subset(
    input_file: [str, Path],
    layer_key: str,
    var_column: str,
    files_to_keep: list,
    slot_map: dict,
    **kwargs,
):
    in_layer = params.get(layer_key, 'X')
    assert in_layer is not None, f'Please specify a layer key in the config under {layer_key}'
    out_layer = f'layers/{layer_key}'
    
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X=in_layer,
        obs='obs',
        var='raw/var' if 'raw/' in in_layer else 'var',
        uns='uns',
        backed=True,
        dask=True,
        **kwargs
    )
    
    logging.info('Subset highly variable genes...')
    subsetted = False
    if var_column not in adata.var:
        adata.var['highly_variable'] = True
    else:
        adata.var['highly_variable'] = adata.var[var_column]
    if save_subset:
        adata, subsetted = subset_hvg(
            adata,
            var_column='highly_variable',
            to_memory=False,
            compute_dask=False,
        )
        dask_compute(adata, layers=[layer_key], verbose=True)
    
    # determine output
    if subsetted and save_subset:
        files_to_keep.append(out_layer)
    else:
        slot_map |= {out_layer: in_layer}
    
    return adata, files_to_keep, slot_map


files_to_keep = ['obs', 'var']
slot_map = {}

adata_norm, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    layer_key='norm_counts',
    var_column=var_mask,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
    obsm='obsm',
    obsp='obsp',
    varm='varm',
    varp='varp',
)
adata_raw, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    layer_key='raw_counts',
    var_column=var_mask,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
)

assert adata_norm.var_names.equals(adata_raw.var_names), f'\nnorm:\n{adata_norm.var}\nraw:\n{adata_raw.var}'
assert adata_norm.n_obs == adata_raw.n_obs, f'\nnorm:\n{adata_norm.obs}\nraw:\n{adata_raw.obs}'

if input_file.endswith('.h5ad'):
    adata = read_anndata(
        input_file,
        obs='obs',
        obsm='obsm',
        obsp='obsp',
        uns='uns',
    )
    adata = AnnData(
        X=adata_norm.X,
        obs=adata.obs,
        var=adata_norm.var,
        obsm=adata.obsm,
        obsp=adata.obsp,
        uns=adata.uns,
        layers={
            'norm_counts': adata_norm.X,
            'raw_counts': adata_raw.X,
        },
    )
    files_to_keep.append('X')
elif input_file.endswith('.zarr'):
    if save_subset:
        files_to_keep.extend(['varm', 'varp'])
    adata = AnnData(
        obs=adata_norm.obs,
        obsm=adata_norm.obsm,
        obsp=adata_norm.obsp,
        var=adata_norm.var,
        layers={
            'norm_counts': adata_norm.X,
            'raw_counts': adata_raw.X,
        },
        uns=adata_norm.uns,
    )
else:
    raise ValueError(f'Invalid input file {input_file}')

# preprocess if missing
if recompute_pca or 'X_pca' not in adata.obsm:
    dask_compute(adata, layers=['norm_counts'], verbose=True)
    logging.info('Compute PCA...')
    import scanpy
    scanpy.pp.pca(adata, layer='norm_counts', mask_var='highly_variable')
    files_to_keep.extend(['obsm', 'uns'])

try:
    assert not recompute_pca, 'PCA was recomputed'
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
except AssertionError as e:
    logging.info(f'Compute neighbors due to Assertion Error: {e}...')
    sc.pp.neighbors(adata, use_rep='X_pca')
    files_to_keep.extend(['obsp', 'uns'])

# fix batch covariate
for batch_key in batch_keys:
    assert batch_key in adata.obs, f'Batch key {batch_key} is missing'
    adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')

# fix labels
for label_key in label_keys:
    if label_key is None or label_key == 'None':
        continue
    assert label_key in adata.obs, f'Label key {label_key} is missing'
    adata.obs[label_key] = adata.obs[label_key].astype(str).astype('category')


logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
slot_map |= {'X': 'layers/norm_counts'}
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
    in_dir_map={
        'layers/norm_counts': output_file
    },
)