from .utils import select_neighbors, rename_categories


def clisi(adata, output_type, meta):
    import scib

    return scib.me.clisi_graph(
        adata,
        batch_key=meta['batch'],
        label_key=meta['label'],
        type_=output_type,
        subsample=0.5 * 100,
        scale=True,
    )


def clisi_y(adata, output_type, meta):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, meta['label'])

    return scib_metrics.clisi_knn(
        X=adata.obsp['distances'],
        labels=labels
    ).mean()


def ilisi(adata, output_type, meta):
    import scib

    return scib.me.ilisi_graph(
        adata,
        batch_key=meta['batch'],
        type_=output_type,
        subsample=0.5 * 100,
        scale=True,
    )


def ilisi_y(adata, output_type, meta):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, meta['batch'])

    return scib_metrics.ilisi_knn(
        X=adata.obsp['distances'],
        batches=batches
    ).mean()
