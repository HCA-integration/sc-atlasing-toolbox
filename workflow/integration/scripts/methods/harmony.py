import logging
logging.basicConfig(level=logging.INFO)
import torch
from harmony import harmonize

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_adata}...')
adata = read_anndata(input_adata)

# prepare for integration
del adata.X

# run method
logging.info('Run harmony...')
adata.obsm["X_emb"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=wildcards.batch)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)