import logging
logging.basicConfig(level=logging.INFO)
import torch
from harmony import harmonize

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    obs='obs',
    var='var',
    obsm='obsm',
    uns='uns'
)

assert 'X_pca' in adata.obsm.keys(), 'PCA is missing'

# run method
logging.info('Run harmony...')
adata.obsm["X_emb"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=wildcards.batch)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])