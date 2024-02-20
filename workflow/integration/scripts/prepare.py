from anndata import AnnData
import scanpy as sc
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.processing import assert_neighbors


def read_layer(
    input_file: [str, Path],
    layer: str,
    **kwargs,
) -> (AnnData, str):
    """
    Read anndata with correct count slot
    :param input_file: input file
    :param layer: slot in file store to read e.g. 'X', 'raw/X', 'layers/counts'
    :return: anndata object, var_key, layer
    """
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X=layer,
        var='raw/var' if 'raw/' in layer else 'var',
        **kwargs,
    )
    return adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params


def read_and_subset(
    input_file: [str, Path],
    layer_key: str,
    files_to_keep: list,
    slot_map: dict,
    hvgs: list = None,
):
    in_layer = params.get(layer_key, 'X')
    out_layer = f'layers/{layer_key}'
    
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X=in_layer,
        var='raw/var' if 'raw/' in in_layer else 'var',
        varm='varm',
        varp='varp',
        backed=True,
    )
    
    logging.info('Subset highly variable genes...')
    adata, subsetted = subset_hvg(adata, hvgs=hvgs)
    
    # determine output
    if subsetted:
        files_to_keep.extend(['X', out_layer, 'var', 'varm', 'varp'])
    else:
        slot_map |= {out_layer: in_layer}
    
    return adata, files_to_keep, slot_map


files_to_keep = []
slot_map = {}

adata_norm, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    layer_key='norm_counts',
    files_to_keep=files_to_keep,
    slot_map=slot_map,
)
adata_raw, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    layer_key='raw_counts',
    files_to_keep=files_to_keep,
    slot_map=slot_map,
    hvgs=adata_norm.var_names,
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
elif input_file.endswith('.zarr'):
    adata = AnnData(
        obs=adata_norm.obs,
        var=adata_norm.var,
        layers={
            'norm_counts': adata_norm.X,
            'raw_counts': adata_raw.X,
        },
    )
else:
    raise ValueError(f'Invalid input file {input_file}')

# preprocess if missing
if 'X_pca' not in adata.obsm:
    logging.info('Compute PCA...')
    adata.X = adata.layers['norm_counts']
    sc.pp.pca(adata)
    del adata.X
    files_to_keep.extend(['obsm', 'varm', 'uns'])

try:
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
except AssertionError:
    logging.info('Compute neighbors...')
    sc.pp.neighbors(adata)
    files_to_keep.extend(['obsp', 'uns'])

logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
)