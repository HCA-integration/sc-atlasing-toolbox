from pathlib import Path
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, remove_slots
from utils.io import read_anndata, write_zarr_linked
from utils.processing import assert_neighbors
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
var_mask = wildcards.var_mask

adata = read_anndata(
    input_file,
    X='layers/norm_counts',
    obs='obs',
    var='var',
    obsp='obsp',
    obsm='obsm',
    varm='varm',
    uns='uns',
    dask=True,
    backed=True,
)

# subset features
adata, _ = subset_hvg(adata, var_column=var_mask)

# prepare output adata
files_to_keep = ['obsm', 'var', 'uns']

if 'X_pca' not in adata.obsm:
    logging.info('Compute PCA...')
    sc.pp.pca(adata)
    files_to_keep.extend(['varm'])
adata.obsm['X_emb'] = adata.obsm['X_pca']

logging.info(adata.__str__())
logging.info(adata.uns.keys())
try:
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
except AssertionError:
    logging.info('Compute neighbors...')
    sc.pp.neighbors(adata)
    print(adata.uns['neighbors'])
    files_to_keep.extend(['obsp', 'uns'])

adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# write file
logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map={'X': 'layers/norm_counts'},
)