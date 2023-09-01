"""
Normalisation
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse
import scanpy as sc

from utils.io import read_anndata, link_zarr


input_file = snakemake.input[0]
output_file = snakemake.output[0]

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)
adata.X = sparse.csr_matrix(adata.X)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

# adata.layers['counts'] = adata.X.copy()
# select counts layer
logging.info('Select layer...')
layer = snakemake.params['raw_counts']
layer = 'X' if layer is None else layer
adata.X = adata.X if layer == 'X' or layer is None else adata.layers[layer]

logging.info('normalize_total...')
sc.pp.normalize_total(adata)
logging.info('log-transform...')
sc.pp.log1p(adata)
adata.X = sparse.csr_matrix(adata.X)

# adata.layers['normcounts'] = adata.X

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = 'default'
adata.uns['preprocessing']['log-transformed'] = True

logging.info(f'Write to {output_file}...')
del adata.raw
del adata.layers
adata.write_zarr(output_file)

# if input_file.endswith('.zarr'):
#     input_files = [f.name for f in Path(input_file).iterdir()]
#     files_to_keep = [f for f in input_files if f not in ['X', 'uns']]
#     link_zarr(
#         in_dir=input_file,
#         out_dir=output_file,
#         file_names=files_to_keep,
#         overwrite=True,
#     )
