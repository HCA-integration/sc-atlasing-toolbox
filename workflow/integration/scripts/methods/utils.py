import anndata as ad
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix, issparse


# TODO: put in common location
def read_anndata(file):
    return ad.read_zarr(file) if file.endswith('.zarr') else sc.read(file)


def select_layer(adata, layer, force_dense=False, force_sparse=False, dtype='float64'):
    # select matrix
    # matrix = adata.X  if layer == 'X' or layer is None else adata.layers[layer]
    if layer == 'X' or layer is None:
        matrix = adata.X
    elif layer in adata.layers:
        matrix = adata.layers[layer]
    elif layer in ['raw', 'counts']:
        try:
            assert not isinstance(adata.raw, type(None))
            matrix = adata.raw[:, adata.var_names].X
        except AssertionError as e:
            raise ValueError(f'Cannot find layer "{layer}" and no counts in adata.raw') from e
    else:
        raise ValueError(f'Invalid layer {layer}')

    if force_dense and force_sparse:
        raise ValueError('force_dense and force_sparse cannot both be True')

    if force_dense:
        matrix = np.asarray(matrix.todense()) if issparse(matrix) else matrix
        return matrix.astype(dtype)

    if force_sparse:
        return matrix if issparse(matrix) else csr_matrix(matrix, dtype=dtype)

    return matrix


def ensure_sparse(adata):
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)


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


def add_metadata(adata, wildcards, params):
    """
    Add integration metatdata to integratd output
    :param adata:
    :param wildcards:
    :param params:
    :return:
    """
    # TODO: transfer parameters from .uns['preprocessing']

    adata.uns['dataset'] = wildcards.dataset

    if 'methods' in adata.uns:
        adata.uns['methods'].append(wildcards.method)
    else:
        adata.uns['methods'] = [wildcards.method]

    adata.uns['integration'] = {
        'method': wildcards.method,
        'label_key': wildcards.label,
        'batch_key': wildcards.batch,
        'output_type': params['output_type'],
        'hyperparams': params['hyperparams'],
    }
