######### UMAP and embedding plots #########

use rule umap from preprocessing as integration_compute_umap with:
    input:
        anndata=rules.integration_postprocess.output.zarr,
        rep=lambda wildcards: get_input_file(config, wildcards, module_name)
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
    params:
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule plot_umap from preprocessing as integration_plot_umap with:
    input:
        anndata=rules.integration_compute_umap.output.zarr,
    output:
        plot=image_dir / 'umap' / f'{paramspace.wildcard_pattern}.png',
        additional_plots=directory(image_dir / 'umap' / paramspace.wildcard_pattern),
    params:
        color=lambda w: [
            get_params(w,parameters,'label'),
            get_params(w,parameters,'batch'),
            *get_for_dataset(
                config=config,
                dataset=w.dataset,
                query=[module_name, 'umap_colors'],
                default=[]
            ),
        ],
        ncols=1,
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
        outlier_factor=10,
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),


rule plots_all:
    input:
        rules.benchmark_all.input,
        expand(rules.integration_plot_umap.output,zip,**parameters[wildcard_names].to_dict('list')),
