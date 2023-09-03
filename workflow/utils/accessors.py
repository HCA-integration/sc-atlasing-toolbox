import numpy as np
from scipy.sparse import csr_matrix, issparse


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


def select_neighbors(adata, output_type):
    neighbors_key = f'neighbors_{output_type}'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    adata.obsp['connectivities'] = adata.obsp[adata.uns[neighbors_key]['connectivities_key']]
    adata.obsp['distances'] = adata.obsp[adata.uns[neighbors_key]['distances_key']]
    return adata
