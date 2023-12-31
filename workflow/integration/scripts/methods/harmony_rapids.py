import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import rapids_singlecell as rsc
import cupy as cp

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', var='var', obsm='obsm', uns='uns')

# run method
logging.info('Run rapids_singlecell harmony...')
rsc.pp.harmony_integrate(
    adata,
    key=wildcards.batch,
    basis='X_pca',
    adjusted_basis='X_emb'
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])