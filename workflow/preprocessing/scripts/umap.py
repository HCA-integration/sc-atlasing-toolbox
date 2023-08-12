"""
Compute UMAP
"""

import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc

from utils.io import read_anndata

input_file = snakemake.input[0]
input_rep = snakemake.input.rep
output_file = snakemake.output[0]
params = {k: v for k, v in snakemake.params.items()}

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write_zarr(output_file)
    exit(0)

# check if representation is available
neighbors_key = params.get('neighbors_key', 'neighbors')
assert neighbors_key in adata.uns
assert adata.uns[neighbors_key]['connectivities_key'] in adata.obsp
assert adata.uns[neighbors_key]['distances_key'] in adata.obsp

# check if representation is available in current anndata, otherwise read from file
use_rep = adata.uns[neighbors_key]['params'].get('use_rep', None)
if use_rep not in adata.obsm:
    logging.info(f'Read {input_rep}...')
    adata_rep = read_anndata(input_rep)
    if use_rep == 'X' or use_rep is None:
        adata.obsm[use_rep] = adata_rep.X
    else:
        adata.obsm[use_rep] = adata_rep.obsm[use_rep]

# compute UMAP
try:
    logging.info('Compute UMAP...')
    sc.tl.umap(adata, method='rapids', **params)
except Exception as e:
    print('sc.tl.umap: Rapids failed, defaulting to UMAP implementation')
    print(e)
    sc.tl.umap(adata, **params)

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)