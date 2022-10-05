import scanpy as sc
import pandas as pd


def prepare_unintegrated(adata, n_hvg=200):
    # TODO: preprocessing according to what is in .uns
    adata_orig = adata.copy()
    sc.pp.highly_variable_genes(adata_orig, n_top_genes=n_hvg)
    sc.pp.pca(adata_orig)
    sc.pp.neighbors(adata_orig)
    return adata_orig


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
