import anndata as ad
import numpy as np


def assert_pca(adata: ad.AnnData):
    assert 'X_pca' in adata.obsm, f'PCA layer missing from obsm:\n{adata}'
    assert 'PCs' in adata.varm, f'PCS missing from varm:\n{adata}'
    assert 'pca' in adata.uns, f'pca missing from uns:\n{adata}'
    assert 'variance' in adata.uns['pca'], f'variance missing from uns[pca]:\n{adata.uns}'
    assert 'variance_ratio' in adata.uns['pca'], f'variance_ratio missing from uns[pca]:\n{adata.uns}'


def assert_neighbors(
    adata,
    neighbors_key='neighbors',
    conn_key='connectivities',
    dist_key='distances',
    check_params=True,
    check_n_neighbors=False,
):
    assert neighbors_key in adata.uns, f'neighbors key "{neighbors_key}" not on .uns'
    assert 'connectivities_key' in adata.uns[neighbors_key]
    assert 'distances_key' in adata.uns[neighbors_key]
    assert adata.uns[neighbors_key]['connectivities_key'] == conn_key, f'"{conn_key} is not saved as conectivities_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert adata.uns[neighbors_key]['distances_key'] == dist_key, f'"{dist_key} is not saved as distances_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert conn_key in adata.obsp, f'"{conn_key}" not in .obsp {adata.obsp.keys()}'
    assert dist_key in adata.obsp, f'"{dist_key}" not in .obsp {adata.obsp.keys()}'
    if check_params:
        assert 'params' in adata.uns[neighbors_key]
        assert 'use_rep' in adata.uns[neighbors_key]['params']
    
    # check that all cells have the same number of neighbors
    n_neighbors = np.unique(adata.obsp['distances'].nonzero()[0], return_counts=True)[1]
    n_neighbors = np.unique(n_neighbors)
    try:
        assert len(n_neighbors) == 1, f'Cells do not have the same number of neighbors {n_neighbors}'
    except AssertionError as e:
        if check_n_neighbors:
            raise e
        print(f'Warning: {e}')