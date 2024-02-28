import scanpy as sc
import pandas as pd


def write_metrics(filename, scores, output_types, **kwargs):
    """
    Write metrics for output type specific scores
    :param filename: file to write to
    :param scores: list of scores, order must match that of output_types
    :param output_types: list of output types, order must match that of scores
    :param kwargs: additional information to add to output
    """
    meta_names = list(kwargs.keys())
    meta_values = [kwargs[col] for col in meta_names]

    records = [
        (*meta_values, output_type, score) 
        for score, output_type in zip(scores, output_types)
    ]
    df = pd.DataFrame.from_records(records, columns=meta_names + ['output_type', 'score'])
    df.to_csv(filename, sep='\t', index=False)


def rename_categories(adata, obs_col):
    s = adata.obs[obs_col]
    s = s.cat.rename_categories({i for i, _ in enumerate(s.cat.categories)})
    return s.to_numpy()


def select_neighbors(adata, output_type):
    # neighbors_key = f'neighbors_{output_type}'
    neighbors_key = 'neighbors'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    
    connectivities_key = adata.uns[neighbors_key]['connectivities_key']
    assert connectivities_key in adata.obsp, f'Connectivities key "{connectivities_key}" missing from adata.obsp {adata}'
    adata.obsp['connectivities'] = adata.obsp[connectivities_key]
    
    distances_key = adata.uns[neighbors_key]['distances_key']
    assert distances_key in adata.obsp, f'Distances key "{distances_key}" missing from adata.obsp {adata}'
    adata.obsp['distances'] = adata.obsp[distances_key]
    
    return adata


# TODO: include in scib package
# TODO: clustering over multiple resolutions in separate rule in parallel
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
