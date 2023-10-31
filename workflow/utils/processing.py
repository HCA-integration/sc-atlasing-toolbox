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


def assert_neighbors(adata, neighbors_key='neighbors', conn_key='connectivities', dist_key='distances', check_params=True):
    assert neighbors_key in adata.uns, f'neighbors key "{neighbors_key}" not on .uns'
    assert adata.uns[neighbors_key]['connectivities_key'] == conn_key, f'"{conn_key} is not saved as conectivities_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert adata.uns[neighbors_key]['distances_key'] == dist_key, f'"{dist_key} is not saved as distances_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert conn_key in adata.obsp, f'"{conn_key}" not in .obsp {adata.obsp.keys()}'
    assert dist_key in adata.obsp, f'"{dist_key}" not in .obsp {adata.obsp.keys()}'
    if check_params:
        assert 'params' in adata.uns[neighbors_key]
        assert 'use_rep' in adata.uns[neighbors_key]['params']


def compute_neighbors(adata, output_type=None, force=False, **kwargs):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :param force: force re-computation of kNN graph
    :param kwargs: additional arguments for sc.pp.neighbors
    :return: anndata with kNN graph based on output type
    """
    neighbors_key = 'neighbors'
    conn_key = 'connectivities'
    dist_key = 'distances'
    
    if not output_type:
        output_type = 'full'
    else:
        neighbors_key = '_'.join([neighbors_key, output_type])
        conn_key = '_'.join([conn_key, output_type])
        dist_key = '_'.join([dist_key, output_type])

    try:
        logging.info(adata.__str__())
        logging.info(adata.uns.get(neighbors_key))
        assert not force
        assert_neighbors(adata, neighbors_key=neighbors_key, conn_key=conn_key, dist_key=dist_key)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError:
        logging.info(adata.__str__())
        logging.info(adata.uns.get(neighbors_key))
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        adata_knn = adata
        if 'params' not in adata_knn.uns['neighbors']:
            adata_knn.uns['neighbors']['params'] = {}
        adata_knn.uns['neighbors']['params'] |= dict(use_rep='X_pca')
        assert_neighbors(adata_knn)

    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm
        adata_knn = adata.copy()
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata_knn, **kwargs)
        assert_neighbors(adata_knn)

    elif output_type == 'full':
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
        adata_knn = adata.copy()
        sc.pp.pca(adata_knn, use_highly_variable=True)
        kwargs |= dict(use_rep='X_pca')
        sc.pp.neighbors(adata_knn, **kwargs)
        assert_neighbors(adata_knn)
        # save PCA
        adata.obsm['X_pca'] = adata_knn.obsm['X_pca']
        adata.varm['PCs'] = adata_knn.varm['PCs']
        adata.uns['pca'] = {
            'variance_ratio':  adata_knn.uns['pca']['variance_ratio'],
            'variance':  adata_knn.uns['pca']['variance'],
        }

    else:
        raise ValueError(f'Invalid output type {output_type}')

    # save kNN graph in dedicated slots
    adata.obsp[conn_key] = adata_knn.obsp['connectivities']
    adata.obsp[dist_key] = adata_knn.obsp['distances']
    adata.uns[neighbors_key] = {
        'connectivities_key': conn_key,
        'distances_key': dist_key,
        'params': adata_knn.uns['neighbors']['params'],
    }
    del adata_knn
    assert_neighbors(adata, neighbors_key, conn_key, dist_key)
