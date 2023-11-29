use rule plot_embedding from preprocessing as preprocessing_plot_embedding with:
    input:
        anndata=rules.pca.output.zarr,
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'pca.png',
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'colors']),
        basis='X_pca',
        ncols=1,


use rule plot_umap from preprocessing as preprocessing_plot_umap with:
    input:
        anndata=rules.umap.output.zarr,
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'umap.png',
        # additional_plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap'),
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'colors']),
        ncols=1,
