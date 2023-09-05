use rule umap from preprocessing as integration_per_lineage_compute_umap with:
    input:
        anndata=rules.integration_per_lineage_postprocess.output.zarr,
        rep=rules.integration_per_lineage_run.input.zarr,
    output:
        zarr=directory(out_dir / 'per_lineage' / 'umap' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'umap.zarr'),
    params:
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule plot_umap from preprocessing as integration_per_lineage_plot_umap with:
    input:
        anndata=rules.integration_per_lineage_compute_umap.output.zarr
    output:
        plot=image_dir / 'umap' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}.png',
        additional_plots=directory(image_dir / 'umap' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}'),
    params:
        color=lambda wildcards: [
            get_params(wildcards,parameters,'label'),
            get_params(wildcards,parameters,'batch'),
            *get_for_dataset(
                config=config,
                dataset=wildcards.dataset,
                query = [module_name, 'umap_colors'],
                default=[]
            ),
        ],
        ncols=1,
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
        outlier_factor=10,
    retries: 2
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


rule umaps_collect:
    input:
        unpack(lambda w: collect_lineages(w, rules.integration_per_lineage_plot_umap.output))
    output:
        touch(image_dir / 'umap' / 'per_lineage' / paramspace.wildcard_pattern / 'umap.done')


rule umaps_all:
    input:
        expand(
            rules.umaps_collect.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        ),