"""
Compute UMAP
"""
from pathlib import Path
import numpy as np
import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')

from utils.io import read_anndata, link_zarr_partial, check_slot_exists
from utils.misc import ensure_dense


def check_and_update_neighbors_info(adata, neighbors_key):
    # check if representation is available
    if neighbors_key not in adata.uns:
        adata.uns[neighbors_key] = {
            'connectivities_key': 'connectivities',
            'distances_key': 'distances',
            'params': {
                'use_rep': 'X_pca',
                'method': None,
            },
        }
    adata.uns[neighbors_key]['params'] = adata.uns[neighbors_key].get('params', {'use_rep': 'X_pca'})
    
    try:
        assert adata.uns[neighbors_key]['connectivities_key'] in adata.obsp
        assert adata.uns[neighbors_key]['distances_key'] in adata.obsp
    except AssertionError as e:
        logging.info(f'neighbors_key: {neighbors_key}')
        logging.info(adata.__str__())
        raise e

    # determine n_neighbors if missing
    if 'n_neighbors' not in adata.uns[neighbors_key]['params']:
        distances_key = adata.uns[neighbors_key]['distances_key']
        n_neighbors = np.unique(adata.obsp[distances_key].nonzero()[0], return_counts=True)[1]
        assert len(np.unique(n_neighbors)) == 1, f'Number of neighbors is not consistent!\nn_neighbors: {np.unique(n_neighbors)}'
        adata.uns[neighbors_key]['params']['n_neighbors'] = int(n_neighbors[0])

    # get correct cell-level representation for UMAP
    use_rep = adata.uns[neighbors_key]['params'].get('use_rep')
    if use_rep not in adata.obsm:
        logging.info(f'Read representation file {input_rep} for use_rep="{use_rep}"...')
        if use_rep == 'X':
            adata = ad.AnnData(
                X=read_anndata(input_rep, X='X').X,
                obs=adata.obs,
                obsp=adata.obsp,
                uns=adata.uns,
            )
            ensure_dense(adata)
        elif check_slot_exists(input_rep, f'obsm/{use_rep}'):
            adata_rep = read_anndata(input_rep, X=f'obsm/{use_rep}')
            ensure_dense(adata_rep)
            adata.obsm[use_rep] = adata_rep.X
        elif use_rep == 'X_pca':
            adata = ad.AnnData(
                X=read_anndata(input_rep, X='X').X,
                obs=adata.obs,
                obsp=adata.obsp,
                uns=adata.uns,
            )
            sc.pp.pca(adata, n_comps=50)
        else:
            raise ValueError(f'Cannot find representation {use_rep} from {input_rep}')
        
    return adata


input_file = snakemake.input[0]
input_rep = snakemake.input.rep
output_file = snakemake.output[0]
params = dict(snakemake.params.items())

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsp='obsp', uns='uns')

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.obsm['X_umap'] = np.zeros((0, 2))
    adata.write_zarr(output_file)
    exit(0)

neighbors_key = params.get('neighbors_key', 'neighbors')
adata = check_and_update_neighbors_info(adata, neighbors_key)

logging.info('Compute UMAP...')
sc.tl.umap(adata, **params)

logging.info(f'Write to {output_file}...')
del adata.X
adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])