import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import harmonypy as hm

from utils import add_metadata, read_anndata, process, select_layer


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

logging.info(f'Read {input_adata}...')
adata_raw = read_anndata(input_adata)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# prepare for integration
adata = adata_raw.copy()
adata.X = select_layer(adata, params['norm_counts'], force_dense=True)

# run method
logging.info('Run harmonypy...')
harmony_out = hm.run_harmony(
    data_mat=adata.obsm['X_pca'],
    meta_data=adata.obs,
    vars_use=[wildcards.batch]
)
adata.obsm['X_emb'] = harmony_out.Z_corr.T

# prepare output adata
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)