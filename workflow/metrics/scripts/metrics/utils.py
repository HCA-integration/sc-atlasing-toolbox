import numpy as np
from scipy import sparse
import scanpy as sc
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO)


def get_from_adata(adata):
    output_type = adata.uns['integration']['output_type']
    output_types = [output_type] if isinstance(output_type, str) else output_type

    return {
        'label': adata.uns['integration']['label_key'],
        'batch': adata.uns['integration']['batch_key'],
        'output_types': output_types
    }


def write_metrics(filename, scores, output_types, lineages, **kwargs):
    """
    Write metrics for output type specific scores
    :param filename: file to write to
    :param scores: list of scores, order must match that of output_types
    :param output_types: list of output types, order must match that of scores
    :param kwargs: additional information to add to output
    """

    meta_names = [key for key in kwargs.keys()]
    meta_values = [kwargs[col] for col in meta_names]

    records = [
        (*meta_values, output_type, score, lineage) 
        for score, output_type, lineage in zip(scores, output_types, lineages)
    ]
    df = pd.DataFrame.from_records(records, columns=meta_names + ['output_type', 'score', 'lineage'])
    df.to_csv(filename, sep='\t', index=False)


def compute_neighbors(adata, output_type):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :return: anndata with kNN graph based on output type
    """
    neighbor_key = f'neighbors_{output_type}'
    conn_key = f'connectivities_{output_type}'
    dist_key = f'distances_{output_type}'

    if neighbor_key in adata.uns and conn_key in adata.obsp and dist_key in adata.obsp:
        logging.info(f'Using pre-computed {output_type} kNN graph')
        return

    if output_type == 'knn':
        adata_knn = adata
        assert 'connectivities' in adata.obsp
        assert 'distances' in adata.obsp

    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm
        adata_knn = adata.copy()
        sc.pp.neighbors(adata_knn, use_rep='X_emb')

    elif output_type == 'full':
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
        adata_knn = adata.copy()
        sc.pp.pca(adata_knn, use_highly_variable=True)
        sc.pp.neighbors(adata_knn, use_rep='X_pca')

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
    adata.uns[neighbor_key] = {
        'connectivities_key': conn_key,
        'distances_key': dist_key,
    }

    adata.obsp[conn_key] = adata_knn.obsp['connectivities']
    adata.obsp[dist_key] = adata_knn.obsp['distances']

    del adata_knn


def select_neighbors(adata, output_type):
    neighbor_key = f'neighbors_{output_type}'
    adata.obsp['connectivities'] = adata.obsp[adata.uns[neighbor_key]['connectivities_key']]
    adata.obsp['distances'] = adata.obsp[adata.uns[neighbor_key]['distances_key']]
    return adata


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
                f'{prefix}{group}': adata[adata.obs[group_key] == group]
                for group in adata.obs[group_key].unique()
            }
        )
    mudata.uns = adata.uns
    return mudata


# TODO: include in scib package
def cluster_optimal(
        adata,
        label_key,
        cluster_key,
        cluster_function=None,
        metric=None,
        resolutions=None,
        use_rep=None,
        inplace=True,
        verbose=False,
        **kwargs
):
    if cluster_function is None:
        cluster_function = sc.tl.leiden

    if metric is None:
        import scib
        metric = scib.me.nmi

    if resolutions is None:
        n = 10
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []

    if use_rep is None:
        try:
            adata.uns['neighbors']
        except KeyError:
            raise RuntimeError('Neighbours must be computed when setting use_rep to None')
    else:
        print(f'Compute neighbors on rep {use_rep}')
        sc.pp.neighbors(adata, use_rep=use_rep)

    for res in resolutions:
        if verbose:
            print(f'resolution: {res}')
        cluster_function(adata, resolution=res, key_added=cluster_key, **kwargs)
        score = metric(adata, label_key, cluster_key)
        score_all.append(score)
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]
        if verbose:
            print(f'optimal score: {score_max}')

    score_all = pd.DataFrame(zip(resolutions, score_all), columns=('resolution', 'score'))

    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering


def rename_categories(adata, obs_col):
    s = adata.obs[obs_col]
    s = s.cat.rename_categories({i for i, _ in enumerate(s.cat.categories)})
    return s.to_numpy()
