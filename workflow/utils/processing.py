import numpy as np
from scipy import sparse
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from .misc import ensure_sparse


def process(adata, output_type, adata_raw=None):
    """
    Process data based on output type.
    If more than one output type is given, use the most processed output type: knn > embed > full
    :param adata: integrated anndata object
    :param adata_raw: anndata object used as input for integration
    :param output_type: string or list of output type
    :return: integrated anndata object with unintegrated anndata in .raw
    """
    if isinstance(output_type, str):
        output_type = [output_type]

    # # save unintegrated count layers
    # adata.layers['counts'] = adata_raw.layers['counts']
    # adata.layers['normcounts'] = adata_raw.layers['normcounts']

    # # save kNN graph of unintegrated object
    # adata.obsp['connectivities_uni'] = adata_raw.obsp['connectivities']
    # adata.obsp['distances_uni'] = adata_raw.obsp['distances']

    # ensure matrix is sparse
    ensure_sparse(adata)

    # remove pre-existing PCA
    if 'X_pca' in adata.obsm:
        del adata.obsm['X_pca']

    if 'full' in output_type:
        pass

    elif 'embed' in output_type:
        assert 'X_emb' in adata.obsm

    elif 'knn' in output_type:
        assert 'connectivities' in adata.obsp
        assert 'distances' in adata.obsp

    else:
        raise ValueError(f'Invalid output type {output_type}')

    # add unintegrated data
    # adata.X = adata_raw.X
    # adata.raw = adata_raw.copy()

    return adata


def compute_neighbors(adata, output_type):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :return: anndata with kNN graph based on output type
    """
    neighbors_key = f'neighbors_{output_type}'
    conn_key = f'connectivities_{output_type}'
    dist_key = f'distances_{output_type}'

    def assert_neighbors(adata, neighbors_key='neighbors', conn_key='connectivities', dist_key='distances', check_params=True):
        assert neighbors_key in adata.uns, f'neighbors key "{neighbors_key}" not on .uns'
        assert adata.uns[neighbors_key]['connectivities_key'] == conn_key, f'"{conn_key} is not saved as conectivities_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
        assert adata.uns[neighbors_key]['distances_key'] == dist_key, f'"{dist_key} is not saved as distances_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
        assert conn_key in adata.obsp, f'"{conn_key}" not in .obsp {adata.obsp.keys()}'
        assert dist_key in adata.obsp, f'"{dist_key}" not in .obsp {adata.obsp.keys()}'
        if check_params:
            assert 'params' in adata.uns[neighbors_key]
            assert 'use_rep' in adata.uns[neighbors_key]['params']

    try:
        logging.info(adata.__str__())
        logging.info(adata.uns.get(neighbors_key))
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
        sc.pp.neighbors(adata_knn, use_rep='X_emb')  # TODO: use rapids
        assert_neighbors(adata_knn)

    elif output_type == 'full':
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
        adata_knn = adata.copy()
        sc.pp.pca(adata_knn, use_highly_variable=True)
        sc.pp.neighbors(adata_knn, use_rep='X_pca')
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


def anndata_to_mudata(adata, group_key, prefix=''):
    import mudata as mu
    import logging
    logging.basicConfig(level=logging.INFO)

    if isinstance(adata, mu.MuData):
        logging.info('Data is already a MuData object')
        mudata = adata
    elif group_key not in adata.obs.columns:
        logging.info('Data is global AnnData object, use generic group name.')
        mudata = mu.MuData({group_key: adata})
    else:
        logging.info('Data is AnnData object, split by group.')
        mudata = mu.MuData(
            {
                f"{prefix}{group.replace(' ', '_').replace('/', '_')}":
                    adata[adata.obs[group_key] == group]
                for group in adata.obs[group_key].unique()
            }
        )
    mudata.uns = adata.uns
    return mudata
