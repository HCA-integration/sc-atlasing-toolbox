from anndata import AnnData
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked, to_memory
from utils.accessors import subset_hvg


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
    adata, subsetted = subset_hvg(adata)
    
    # determine output
    if subsetted:
        files_to_keep.append(out_layer)
    else:
        slot_map |= {out_layer: in_layer}
    
    return adata, files_to_keep, slot_map


files_to_keep = []
slot_map = {}

adata_norm, files_to_link, slot_map = read_and_subset(
    input_file, 'norm_counts', files_to_keep, slot_map
)
adata_raw, files_to_link, slot_map = read_and_subset(
    input_file, 'raw_counts', files_to_keep, slot_map
)

assert adata_norm.var.equals(adata_raw.var)
assert adata_norm.n_obs == adata_raw.n_obs, f'norm: {adata_norm.n_obs} , raw: {adata_raw.n_obs}'

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

logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
)