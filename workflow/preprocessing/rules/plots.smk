use rule plots from preprocessing as preprocessing_plot_pca with:
    input:
        anndata=rules.pca.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'pca'),
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'colors']),
        basis='X_pca',
        ncols=1,


use rule plots from preprocessing as preprocessing_plot_umap with:
    input:
        anndata=rules.umap.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap'),
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'colors']),
        basis='X_umap',
        ncols=1,
        outlier_factor=100,
