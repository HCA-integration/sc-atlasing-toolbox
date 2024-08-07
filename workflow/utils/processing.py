import traceback
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse
from tqdm import tqdm
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    logging.info('Using rapids_singlecell...')
    USE_GPU = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    USE_GPU = False

from .assertions import assert_neighbors
from .io import to_memory, csr_matrix_int64_indptr


def compute_neighbors(adata, output_type=None, force=False, check_n_neighbors=False, **kwargs):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output typefilter
    :param force: force re-computation of kNN graph
    :param kwargs: additional arguments for sc.pp.neighbors
    :return: anndata with kNN graph based on output type
    """
    
    if not output_type:
        output_type = 'full'

    try:
        logging.info(adata.__str__())
        assert not force, 'force neighbors computation'
        assert_neighbors(adata, check_n_neighbors=check_n_neighbors)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError as e:
        logging.info(e.__str__())
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        if 'params' not in adata.uns['neighbors']:
            adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
        adata.uns['neighbors']['params'] |= dict(use_rep='X')
    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm, 'Embedding key "X_emb" not found'
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata, **kwargs)
    elif output_type == 'full':
        if 'X_emb' not in adata.obsm:
            logging.info('Compute PCA on corrected counts...')
            adata_tmp = adata.copy()
            sc.pp.pca(
                adata_tmp,
                use_highly_variable='highly_variable' in adata_tmp.var.columns
            )
            adata.obsm['X_emb'] = adata_tmp.obsm['X_pca']
            del adata_tmp
        kwargs |= dict(use_rep='X_emb')
        logging.info('Neighbors...')
        sc.pp.neighbors(adata, **kwargs)
    else:
        raise ValueError(f'Invalid output type {output_type}')
    
    assert_neighbors(adata, check_n_neighbors=False)


def _filter_genes(adata, **kwargs):
    import scanpy as sc
    from dask import array as da
    
    if isinstance(adata.X, da.Array):
        adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
    gene_subset, _ = sc.pp.filter_genes(adata.X, **kwargs)
    if isinstance(gene_subset, da.Array):
        gene_subset = gene_subset.compute()
    return gene_subset


def filter_genes(adata, batch_key=None, **kwargs):
    """
    Filter anndata based on .X matrix
    """
    from dask import array as da
    
    # filter cells from batches with too few cells
    if batch_key is not None:
        cells_per_batch = adata.obs[batch_key].value_counts()
        if cells_per_batch.min() < 2:
            adata = adata[adata.obs[batch_key].isin(cells_per_batch[cells_per_batch > 1].index)]
        
    # apply filter
    gene_subset = _filter_genes(adata, **kwargs)
    if any(gene_subset == False):
        logging.info(f'Subset to {sum(gene_subset)}/{adata.n_vars} filtered genes...')
        adata = adata[:, gene_subset]
    
    if adata.is_view:
        adata = adata.copy()    
    
    # convert dask array to csr matrix
    if isinstance(adata.X, da.Array):
        adata.X = adata.X.map_blocks(csr_matrix_int64_indptr).compute()

    return adata


def get_pseudobulks(adata, group_key, agg='sum'):
    from dask import array as da
    
    def aggregate(x, agg):
        if agg == 'sum':
            return x.sum(0)
        elif agg == 'mean':
            return x.mean(0)
        else:
            raise ValueError(f'invalid aggregation method "{agg}"')
    
    def _get_pseudobulk_matrix(adata, group_key, agg):
        X = adata.X
        value_counts = adata.obs[group_key].value_counts()
        
        # filter groups for at least 2 replicates
        value_counts = value_counts[value_counts >= 2]
        groups = value_counts.index
        
        print(f'Aggregate {len(groups)} pseudobulks...')
        
        if isinstance(X, da.Array):
            from tqdm.dask import TqdmCallback
            
            # sort cells by group_key
            df = adata.obs.reset_index(drop=True).query(f'{group_key} in @groups')
            df[group_key] = pd.Categorical(df[group_key], categories=groups, ordered=True)
            sorted_groups = df.sort_values(by=group_key)
            
            # sort dask array by group_key and rechunk by size of groups
            # the result should be a dask chunk per pseudobulk group
            X = X[sorted_groups.index.values].rechunk((tuple(value_counts.values), -1))
            pseudobulks = X.map_blocks(lambda x: aggregate(x, agg))
            
            miniters = max(10, len(pseudobulks.dask) // 100)
            mininterval = min(max(1, len(groups) / 500), 10)
            with TqdmCallback(desc="Aggregate Dask array", miniters=miniters, mininterval=mininterval):
                pseudobulks = pseudobulks.compute()

        elif isinstance(X, (scipy.sparse.spmatrix, np.ndarray)):
            import scanpy as sc

            pbulk_adata = sc.get.aggregate(
                adata,
                by=group_key,
                func=[agg]
            )
            pbulk_adata = pbulk_adata[pbulk_adata.obs[group_key].isin(groups)].copy()
            pseudobulks = pbulk_adata.layers[agg]
            groups = pbulk_adata.obs_names
            # miniters = max(10, len(value_counts) // 100)
            # pseudobulks = []
            # for group in tqdm(value_counts.index, desc='Aggregate groups', miniters=miniters):
            #     row_agg = aggregate(adata[adata.obs[group_key] == group].X, agg)
            #     row_agg = row_agg.A1 if isinstance(row_agg, np.matrix) else row_agg
            #     pseudobulks.append(row_agg)
            # pseudobulks = np.stack(pseudobulks, axis=0)
        else:
            raise ValueError(f'invalid type "{type(x)}"')
        
        return scipy.sparse.csr_matrix(pseudobulks), groups


    pbulks, groups = _get_pseudobulk_matrix(adata, group_key=group_key, agg=agg)

    print(f'Merge metadata...', flush=True)
    obs = pd.DataFrame(groups, columns=[group_key]).merge(
        adata.obs.drop_duplicates(subset=[group_key]),
        on=group_key,
        how='left'
    )
    obs.index = obs[group_key].values

    for col in obs.columns:
        if obs[col].dtype.name == 'category':
            obs[col] = obs[col].astype(str).astype('category')

    return ad.AnnData(
        pbulks,
        obs=obs,
        var=adata.var,
        dtype='float32'
    )
