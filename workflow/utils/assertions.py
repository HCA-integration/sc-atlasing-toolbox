import anndata as ad


def assert_pca(adata: ad.AnnData):
    assert 'X_pca' in adata.obsm, f'PCA layer missing from obsm:\n{adata}'
    assert 'PCs' in adata.varm, f'PCS missing from varm:\n{adata}'
    assert 'pca' in adata.uns, f'pca missing from uns:\n{adata}'
    assert 'variance' in adata.uns['pca'], f'variance missing from uns[pca]:\n{adata.uns}'
    assert 'variance_ratio' in adata.uns['pca'], f'variance_ratio missing from uns[pca]:\n{adata.uns}'
