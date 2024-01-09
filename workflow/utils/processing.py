import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from scipy import sparse
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')


from .assertions import assert_neighbors
from .io import to_memory


def compute_neighbors(adata, output_type=None, force=False, **kwargs):
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
        assert not force
        assert_neighbors(adata)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError:
        logging.info(adata.__str__())
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        if 'params' not in adata.uns['neighbors']:
            adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
        adata.uns['neighbors']['params'] |= dict(use_rep='X')
    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata, **kwargs)
    elif output_type == 'full':
        adata.X = to_memory(adata.X)
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
        kwargs |= dict(use_rep='X')
        sc.pp.pca(adata, use_highly_variable=True)
        sc.pp.neighbors(adata, **kwargs)
    else:
        raise ValueError(f'Invalid output type {output_type}')
    
    assert_neighbors(adata)
