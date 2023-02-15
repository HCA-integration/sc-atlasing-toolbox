from .utils import select_neighbors


def graph_connectivity(adata, output_type, meta):
    import scib

    adata = select_neighbors(adata, output_type)
    return scib.me.graph_connectivity(
        adata,
        label_key=meta['label']
    )
