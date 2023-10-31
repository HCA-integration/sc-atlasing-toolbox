import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import harmonypy as hm

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

adata.write_zarr(output_adata)