import numpy as np
from .utils import select_neighbors, rename_categories


def clisi(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    return scib.me.clisi_graph(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        type_=output_type,
        subsample=0.5 * 100,
        scale=True,
    )


def clisi_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, label_key)
    
    return scib_metrics.clisi_knn(
        X=adata.obsp['distances'],
        labels=labels
    ).mean()


def ilisi(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    return scib.me.ilisi_graph(
        adata,
        batch_key=batch_key,
        type_=output_type,
        subsample=0.5 * 100,
        scale=True,
    )


def ilisi_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, batch_key)

    return scib_metrics.ilisi_knn(
        X=adata.obsp['distances'],
        batches=batches
    ).mean()


def ilisi_re(adata, output_type, batch_key, label_key, **kwargs):
    return lisi_re(adata, output_type, batch_key)


def clisi_re(adata, output_type, batch_key, label_key, **kwargs):
    return lisi_re(adata, output_type, label_key)


def lisi_re(adata, output_type, obs_key):
    import anndata
    import scanpy as sc
    from sklearn.neighbors import NearestNeighbors, sort_graph_by_row_values
    from numba import njit
    
    
    def compute_lisi(
        adata: anndata.AnnData,
        obs_column: str,
        perplexity: float=30,
        use_rep: str=None,
        recompute_knn: bool=True,
    ):
        """
        Code adapted from https://github.com/slowkow/harmonypy/blob/182a5c6e0fc954cb0b4db5074e507adb7a8293f3/harmonypy/lisi.py#L24
        """
        if use_rep is None:
            use_rep = 'X'
            
        labels = adata.obs[obs_column].astype('category')
        labels = labels.map({label: i for i, label in enumerate(labels.cat.categories)})
        labels_categories = np.array(labels.unique())
        labels = np.array(labels)

        # We need at least 3 * n_neighbors to compute the perplexity
        n_neighbors = perplexity * 3
        if recompute_knn:
            print(f'recompute kNN graph with {n_neighbors} neighbors...')
            if use_rep != 'X' and use_rep not in adata.obsm:
                adata.obsm[use_rep] = adata.obsp[use_rep]
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
        
        # check if kNN degree distribution is correct
        inferred_degree = np.unique(np.unique(adata.obsp['distances'].nonzero()[0], return_counts=True)[1])
        assert len(inferred_degree) == 1, f'inferred degrees: {inferred_degree}'
        assert inferred_degree[0] == n_neighbors-1, f'inferred: {inferred_degree[0]}, required degree: {n_neighbors-1}'
        
        # retrieve kNN distances and indices
        X = adata.obsp['distances']
        print('sort graph...')
        X = sort_graph_by_row_values(X)
        knn = NearestNeighbors(n_neighbors=n_neighbors-1, metric="precomputed").fit(X)
        distances, indices = knn.kneighbors(X)

        indices = indices[:,1:]
        distances = distances[:,1:]

        simpson = compute_simpson_matrix(
            distances.T,
            indices.T,
            labels,
            labels_categories,
            perplexity
        )
        return 1 / simpson
    
    @njit
    def compute_simpson_matrix(
        distances: np.ndarray,
        indices: np.ndarray,
        labels: np.ndarray,
        labels_categories: np.ndarray,
        perplexity: float,
        tol: float=1e-5
    ):
        """
        Code adapted from: https://github.com/slowkow/harmonypy/blob/182a5c6e0fc954cb0b4db5074e507adb7a8293f3/harmonypy/lisi.py#L68
        """
        n = distances.shape[1]
        P = np.zeros(distances.shape[0])
        simpson = np.zeros(n)
        logU = np.log(perplexity)
        # Loop through each cell.
        for i in range(n):
            beta = 1
            betamin = -np.inf
            betamax = np.inf
            # Compute Hdiff
            P = np.exp(-distances[:,i] * beta)
            P_sum = np.sum(P)
            if P_sum == 0:
                H = 0
                P = np.zeros(distances.shape[0])
            else:
                H = np.log(P_sum) + beta * np.sum(distances[:,i] * P) / P_sum
                P = P / P_sum
            Hdiff = H - logU
            n_tries = 50
            for t in range(n_tries):
                # Stop when we reach the tolerance
                if abs(Hdiff) < tol:
                    break
                # Update beta
                if Hdiff > 0:
                    betamin = beta
                    if not np.isfinite(betamax):
                        beta *= 2
                    else:
                        beta = (beta + betamax) / 2
                else:
                    betamax = beta
                    if not np.isfinite(betamin):
                        beta /= 2
                    else:
                        beta = (beta + betamin) / 2
                # Compute Hdiff
                P = np.exp(-distances[:,i] * beta)
                P_sum = np.sum(P)
                if P_sum == 0:
                    H = 0
                    P = np.zeros(distances.shape[0])
                else:
                    H = np.log(P_sum) + beta * np.sum(distances[:,i] * P) / P_sum
                    P = P / P_sum
                Hdiff = H - logU
            # distancesefault value
            if H == 0:
                simpson[i] = -1
            # Simpson's index
            for label_category in labels_categories:
                ix = indices[:,i]
                if label_category in labels[ix]:
                    P_sum = np.sum(P[labels[ix] == label_category])
                    simpson[i] += P_sum * P_sum
        
        return simpson
    
    use_rep = 'X_pca'
    if output_type == 'knn':
        use_rep = 'distances'
    elif output_type == 'embed':
        use_rep = 'X_emb'
    
    lisi_scores = compute_lisi(
        adata,
        obs_key,
        perplexity=30,
        use_rep=use_rep,
        recompute_knn=True
    )
    lisi = np.nanmedian(lisi_scores)
    n_labels = adata.obs[obs_key].nunique()
    lisi = (lisi - 1) / (n_labels - 1)
    return lisi
