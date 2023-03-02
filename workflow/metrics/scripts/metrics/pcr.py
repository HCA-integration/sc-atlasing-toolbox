import numpy as np


def pcr(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    adata_raw = adata.raw.to_adata()

    return scib.me.pcr_comparison(
        adata_pre=adata_raw,
        adata_post=adata,
        covariate=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def cell_cycle(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        np.nan

    adata_raw = adata.raw.to_adata()

    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']
        adata_raw.var_names = adata_raw.var['feature_name']

    print(adata.var)
    print(adata.var_names)

    return scib.me.cell_cycle(
        adata_pre=adata_raw,
        adata_post=adata,
        batch_key=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
        organism='human',
    )
