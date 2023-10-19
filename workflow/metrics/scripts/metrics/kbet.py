from .utils import rename_categories, select_neighbors


def kbet_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, batch_key)

    score, _ , _ = scib_metrics.kbet(
        X=adata.obsp['distances'],
        batches=batches,
    )

    return score