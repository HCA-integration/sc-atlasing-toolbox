"""
PCA on highly variable genes
"""

import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse
import scanpy as sc

from utils.io import read_anndata


input_file = snakemake.input[0]
input_counts = snakemake.input.counts
output_file = snakemake.output[0]
scale = snakemake.params['scale']

logging.info(f'Read "{input_file}"...')
adata = read_anndata(input_file)
adata.X = read_anndata(input_counts)[:, adata.var_names].X

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

# scaling
if scale:
    logging.info('Scale data...')
    sc.pp.scale(adata, max_value=10, zero_center=True)
else:
    adata.X = sparse.csr_matrix(adata.X)

logging.info('PCA...')
sc.pp.pca(adata, use_highly_variable=True)

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['scaled'] = scale

# remove counts
del adata.X
del adata.layers

logging.info(f'Write to "{output_file}"...')
adata.write_zarr(output_file)
