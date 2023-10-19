use rule plot_embedding from preprocessing as preprocessing_plot_embedding with:
    input:
        anndata=rules.pca.output.zarr,
    output:
        plot=image_dir / wildcard_pattern / 'pca.png',
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'colors']),
        basis='X_pca',
        ncols=1,


use rule plot_umap from preprocessing as preprocessing_plot_umap with:
    input:
        anndata=rules.umap.output.zarr,
    output:
        plot=image_dir / wildcard_pattern / 'umap.png',
        additional_plots=directory(image_dir / wildcard_pattern / 'umap'),
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'colors']),
        ncols=1,
        outlier_factor=10
