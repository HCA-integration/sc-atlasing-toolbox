import scanpy as sc
from scipy.sparse import csr_matrix


def process(adata, adata_raw, output_type):
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

    if 'full' in output_type:
        sc.pp.pca(adata)
        # TODO: use same parameters as in .uns['preprocessing']

    elif 'embed' in output_type or 'knn' in output_type:
        # remove unintegrated entries for embed and knn
        adata.X = csr_matrix((adata.n_obs, adata.n_vars), dtype='float32')
        if 'X_pca' in adata.obsm:
            del adata.obsm['X_pca']

    else:
        raise ValueError(f'Invalid output type {output_type}')

    # save unintegrated count layers
    adata.layers['counts'] = adata_raw.layers['counts']
    adata.layers['normcounts'] = adata_raw.layers['normcounts']
    adata.raw = adata_raw
    return adata
