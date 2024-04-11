import traceback
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import anndata as ad
from scipy import sparse
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    logging.info('Using rapids_singlecell...')
    USE_GPU = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    USE_GPU = False

from .assertions import assert_neighbors
from .io import to_memory


def compute_neighbors(adata, output_type=None, force=False, check_n_neighbors=False, **kwargs):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :param force: force re-computation of kNN graph
    :param kwargs: additional arguments for sc.pp.neighbors
    :return: anndata with kNN graph based on output type
    """
    
    if not output_type:
        output_type = 'full'

    try:
        logging.info(adata.__str__())
        assert not force, 'force neighbors computation'
        assert_neighbors(adata, check_n_neighbors=check_n_neighbors)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError as e:
        logging.info(e.__str__())
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        if 'params' not in adata.uns['neighbors']:
            adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
        adata.uns['neighbors']['params'] |= dict(use_rep='X')
    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm, 'Embedding key "X_emb" not found'
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata, **kwargs)
    elif output_type == 'full':
        if 'X_emb' not in adata.obsm:
            logging.info('Compute PCA on corrected counts...')
            adata_tmp = adata.copy()
            sc.pp.pca(
                adata_tmp,
                use_highly_variable='highly_variable' in adata_tmp.var.columns
            )
            adata.obsm['X_emb'] = adata_tmp.obsm['X_pca']
            del adata_tmp
        kwargs |= dict(use_rep='X_emb')
        logging.info('Neighbors...')
        sc.pp.neighbors(adata, **kwargs)
    else:
        raise ValueError(f'Invalid output type {output_type}')
    
    assert_neighbors(adata, check_n_neighbors=False)
