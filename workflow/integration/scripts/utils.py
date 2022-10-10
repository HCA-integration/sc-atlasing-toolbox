import scanpy as sc
from scipy.sparse import csr_matrix


def process(adata, adata_raw, output_type):
    if isinstance(output_type, str):
        output_type = [output_type]

    if 'knn' in output_type:
        # remove unintegrated entries
        adata.X = csr_matrix((adata.n_obs, adata.n_vars), dtype='float32')
        if 'X_pca' in adata.obsm:
            del adata.obsm['X_pca']

        # assert integrated entries
        assert 'connectivities' in adata.obsp
        assert 'distances' in adata.obsp

    elif 'embed' in output_type:
        # remove unintegrated entries
        adata.X = csr_matrix((adata.n_obs, adata.n_vars), dtype='float32')
        if 'X_pca' in adata.obsm:
            del adata.obsm['X_pca']

        # process integrated entries
        assert 'X_emb' in adata.obsm
        sc.pp.neighbors(adata, use_rep='X_emb')

    elif 'full' in output_type:
        sc.pp.pca(adata)
        adata.obsm['X_emb'] = adata.obsm['X_pca']
        sc.pp.neighbors(adata, use_rep='X_emb')

    else:
        raise ValueError(f'Invalid output type {output_type}')

    adata.layers['counts'] = adata_raw.layers['counts']
    adata.layers['normcounts'] = adata_raw.layers['normcounts']
    adata.raw = adata_raw
    return adata
