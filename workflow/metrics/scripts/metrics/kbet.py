from .utils import rename_categories, select_neighbors


def kbet_y(adata, output_type, meta):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, meta['batch'])

    score, _ , _ = scib_metrics.kbet(
        X=adata.obsp['distances'],
        batches=batches,
    )

    return score