import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import rapids_singlecell as rsc
import cupy as cp

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer

input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

logging.info(f'Read {input_adata}...')
adata = read_anndata(input_adata)

# prepare for integration
del adata.X

# run method
logging.info('Run rapids_singlecell harmony...')
rsc.tl.harmony_integrate(
    adata,
    key=wildcards.batch,
    basis='X_pca',
    adjusted_basis='X_emb'
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)