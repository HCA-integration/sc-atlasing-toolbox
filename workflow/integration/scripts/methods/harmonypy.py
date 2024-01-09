import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import harmonypy as hm

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

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
logging.info('Run harmonypy...')
harmony_out = hm.run_harmony(
    data_mat=adata.obsm['X_pca'],
    meta_data=adata.obs,
    vars_use=[wildcards.batch]
)
adata.obsm['X_emb'] = harmony_out.Z_corr.T

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)