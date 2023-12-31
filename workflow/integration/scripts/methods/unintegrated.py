from pathlib import Path
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial
from utils_pipeline.accessors import select_layer
from utils_pipeline.processing import assert_neighbors


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers', obsp='obsp', uns='uns')
adata.X = select_layer(adata, params['norm_counts'])

# prepare output adata
files_to_keep = ['obsm', 'uns']

if 'X_pca' not in adata.obsm:
    sc.pp.pca(adata, use_highly_variable=True)
    files_to_keep.extend(['varm'])
adata.obsm['X_emb'] = adata.obsm['X_pca']

logging.info(adata.__str__())
logging.info(adata.uns.keys())
try:
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
    files_to_keep.extend(['obsp', 'uns'])
except AssertionError:
    logging.info('Compute neighbors...')
    sc.pp.neighbors(adata)
    print(adata.uns['neighbors'])

adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# write file
adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=files_to_keep)
