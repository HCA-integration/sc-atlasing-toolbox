import anndata as ad


def add_wildcards(
    adata: ad.AnnData,
    wildcards: dict,
    module_name: str,
    exclude: list = ['dataset', 'file_id']
):
    """
    Add wildcards to adata.uns
    """
    adata.uns['wildcards'] = adata.uns.get('wildcards', {})
    adata.uns['wildcards'] |= {
        f'{module_name}_{k}': v for k, v in wildcards.items() if k not in exclude
    }