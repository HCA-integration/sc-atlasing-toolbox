import warnings
import numpy as np
import anndata as ad
from dask import array as da

from .io import to_memory
from .misc import apply_layers


# deprecated
def select_layer(adata, layer, force_dense=False, force_sparse=False, dtype='float32'):
    from scipy.sparse import csr_matrix, issparse
    from dask.array import Array as DaskArray
    
    # select matrix
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

    if isinstance(matrix, DaskArray):
        return matrix

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


def subset_hvg(
    adata: ad.AnnData,
    to_memory: [str, list] = 'X',
    var_column: str = 'highly_variable'
) -> (ad.AnnData, bool):
    """
    Subset to highly variable genes
    :param adata: anndata object
    :param to_memory: layers to convert to memory
    :return: subsetted anndata object, bool indicating whether subset was performed
    """
    assert var_column in adata.var.columns, f'Column {var_column} not found in adata.var'
    assert adata.var[var_column].dtype == bool, f'Column {var_column} is not boolean'

    if adata.var[var_column].sum() == adata.var.shape[0]:
        warnings.warn('All genes are highly variable, not subsetting')
        subsetted = False
    else:
        subsetted = True
        adata = adata[:, adata.var[var_column]].copy()
    adata = adata_to_memory(adata, layers=to_memory)
    adata = apply_layers(
        adata,
        func=lambda x: x.compute() if isinstance(x, da.Array) else x,
        layers=to_memory,
    )
    return adata, subsetted


def adata_to_memory(adata: ad.AnnData, layers: [str, list] = None) -> ad.AnnData:
    if layers is None:
        layers = ['X', 'raw'] + list(adata.layers.keys())
    elif isinstance(layers, str):
        layers = [layers]
    
    for layer in layers:
        if layer in adata.layers:
            adata.layers[layer] = to_memory(adata.layers[layer])
        elif layer == 'X':
            adata.X = to_memory(adata.X)
        elif layer == 'raw':
            if adata.raw is None:
                continue
            adata_raw = adata.raw.to_adata()
            adata_raw.X = to_memory(adata.raw.X)
            adata.raw = adata_raw
    return adata
