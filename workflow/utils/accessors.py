import warnings
import numpy as np
import anndata as ad
from dask import array as da

from .io import to_memory
from .misc import dask_compute
from .processing import _filter_genes


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
    to_memory: [str, list, bool] = 'X',
    var_column: str = 'highly_variable',
    compute_dask: bool = True,
    add_column: str = 'highly_variable',
    min_cells: int = 0,
) -> (ad.AnnData, bool):
    """
    Subset to highly variable genes
    :param adata: anndata object
    :param to_memory: layers to convert to memory
    :return: subsetted anndata object, bool indicating whether subset was performed
    """
    import sparse
    from tqdm.dask import TqdmCallback
    
    if var_column is None or var_column == 'None':
        var_column = 'mask'
        adata.var[var_column] = True
    assert var_column in adata.var.columns, f'Column {var_column} not found in adata.var'
    assert adata.var[var_column].dtype == bool, f'Column {var_column} is not boolean'
    
    if add_column is not None and add_column != var_column:
        adata.var[add_column] = adata.var[var_column]
    
    if min_cells > 0:
        print(f'Filter to genes that are in at least {min_cells} cells...', flush=True)
        if isinstance(adata.X, da.Array):
            low_count_mask = adata.X.map_blocks(sparse.COO).sum(axis=0) < min_cells
            with TqdmCallback(desc=f"Determine genes with < {min_cells} cells"):
                low_count_mask = low_count_mask.compute().todense()
        else:
            low_count_genes = _filter_genes(adata, min_cells=min_cells)
            low_count_mask = np.isin(adata.var_names, low_count_genes)
        adata.var[var_column] = adata.var[var_column] & ~low_count_mask
    
    if adata.var[var_column].sum() == adata.var.shape[0]:
        warnings.warn('All genes are highly variable, not subsetting')
        subsetted = False
    else:
        subsetted = True
        print(
            f'Subsetting to {adata.var[var_column].sum()} features from {var_column}...',
            flush=True
        )
        adata._inplace_subset_var(adata.var_names[adata.var[var_column]])
    
    if compute_dask:
        adata = dask_compute(adata, layers=to_memory)
    else:
        adata = adata_to_memory(
            adata,
            layers=to_memory,
            verbose=False,
        )
    
    return adata, subsetted


def adata_to_memory(
    adata: ad.AnnData,
    layers: [str, list, bool] = None,
    verbose: bool = False,
) -> ad.AnnData:
    if layers is None or layers is True:
        layers = ['X', 'raw'] + list(adata.layers.keys())
    elif isinstance(layers, str):
        layers = [layers]
    elif layers is False:
        return adata
    
    for layer in layers:
        if verbose:
            print(f'Convert {layer} to memory...', flush=True)
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
        elif verbose:
            print(f'Layer {layer} not found, skipping...', flush=True)
    return adata
