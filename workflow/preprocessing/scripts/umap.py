"""
Compute UMAP
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from anndata.experimental import read_elem
import zarr
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')

from utils.io import read_anndata, link_zarr
from utils.misc import ensure_dense


def check_and_update_neighbors_info(adata, neighbors_key):
    # check if representation is available
    try:
        if neighbors_key not in adata.uns:
            adata.uns[neighbors_key] = {
                'connectivities_key': 'connectivities',
                'distances_key': 'distances',
                'params': {
                    'use_rep': 'X_pca',
                    'method': None,
                },
            }
        assert adata.uns[neighbors_key]['connectivities_key'] in adata.obsp
        assert adata.uns[neighbors_key]['distances_key'] in adata.obsp
    except AssertionError as e:
        logging.info(f'neighbors_key: {neighbors_key}')
        logging.info(adata.__str__())
        raise e

    # check if representation is available in current anndata, otherwise read from file
    use_rep = adata.uns[neighbors_key]['params'].get('use_rep', None)
    if use_rep not in adata.obsm:
        logging.info(f'Read representation file {input_rep}...')
        use_counts = (use_rep == 'X') or (use_rep is None)
        if use_counts:
            adata.X = read_anndata(input_rep, X='X').X
            ensure_dense(adata)
        else:
            adata.obsm[use_rep] = read_anndata(input_rep, obs='obs', obsm='obsm').obsm[use_rep]


input_file = snakemake.input[0]
input_rep = snakemake.input.rep
output_file = snakemake.output[0]
params = dict(snakemake.params.items())

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm', obsp='obsp', uns='uns', var='var')

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.obsm['X_umap'] = np.zeros((0, 2))
    adata.write_zarr(output_file)
    exit(0)

neighbors_key = params.get('neighbors_key', 'neighbors')

# compute UMAP
# for key in neighbors_key:
#     check_and_update_neighbors_info(adata, key)
#     params |= dict(neighbors_key=key)
#     sc.tl.umap(adata, **params)
#     if key != 'neighbors' or len(neighbors_key) > 1:
#         adata.obsm[f'X_umap_{key}'] = adata.obsm['X_umap']

check_and_update_neighbors_info(adata, neighbors_key)
sc.tl.umap(adata, **params)

logging.info(f'Write to {output_file}...')
del adata.X
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obsm', 'uns']]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_keep,
        overwrite=True,
)