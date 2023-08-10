from .utils import select_neighbors, rename_categories


def graph_connectivity(adata, output_type, meta, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    return scib.me.graph_connectivity(
        adata,
        label_key=meta['label']
    )


def graph_connectivity_y(adata, output_type, meta, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, meta['label'])

    return scib_metrics.clisi_knn(
        X=adata.obsp['distances'],
        labels=labels
    )
