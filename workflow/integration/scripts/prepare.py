from anndata import AnnData
from dask import array as da
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import apply_layers
from utils.processing import assert_neighbors, sc


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
var_mask = snakemake.wildcards.var_mask
save_subset = params.get('save_subset', False)


def read_and_subset(
    input_file: [str, Path],
    layer_key: str,
    var_column: str,
    files_to_keep: list,
    slot_map: dict,
):
    in_layer = params.get(layer_key, 'X')
    assert in_layer is not None, f'Please specify a layer key in the config under {layer_key}'
    out_layer = f'layers/{layer_key}'
    
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X=in_layer,
        obsm='obsm',
        var='raw/var' if 'raw/' in in_layer else 'var',
        varm='varm',
        varp='varp',
        uns='uns',
        backed=True,
        dask=True,
    )
    
    logging.info('Subset highly variable genes...')
    adata, subsetted = subset_hvg(
        adata,
        var_column=var_column,
        to_memory=False,
        compute_dask=False,
    )
    
    # determine output
    if subsetted and save_subset:
        files_to_keep.extend([out_layer, 'var', 'varm', 'varp'])
    else:
        slot_map |= {out_layer: in_layer}
    
    return adata, files_to_keep, slot_map


files_to_keep = []
slot_map = {}

adata_norm, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    layer_key='norm_counts',
    var_column=var_mask,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
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
    adata = AnnData(
        obs=adata_norm.obs,
        obsm=adata_norm.obsm,
        var=adata_norm.var,
        layers={
            'norm_counts': adata_norm.X,
            'raw_counts': adata_raw.X,
        },
        uns=adata_norm.uns,
    )
else:
    raise ValueError(f'Invalid input file {input_file}')

apply_layers(
    adata,
    func=lambda x: x.compute() if isinstance(x, da.Array) else x,
    layers=['norm_counts', 'raw_counts'],
    verbose=True,
)

# preprocess if missing
if 'X_pca' not in adata.obsm:
    logging.info('Compute PCA...')
    import scanpy
    scanpy.pp.pca(adata, layer='norm_counts')
    files_to_keep.extend(['obsm', 'uns'])

try:
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
except AssertionError:
    logging.info('Compute neighbors...')
    sc.pp.neighbors(adata)
    files_to_keep.extend(['obsp', 'uns'])

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