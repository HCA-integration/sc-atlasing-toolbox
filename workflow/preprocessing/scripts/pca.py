"""
PCA on highly variable genes
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse
import scanpy as sc

from utils.io import read_anndata, link_zarr


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

logging.info(f'Write to "{output_file}"...')
del adata.raw
del adata.X
del adata.layers
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_link = [f for f in input_files if f not in ['obsm', 'uns', 'varm']]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_link,
        overwrite=True,
    )
