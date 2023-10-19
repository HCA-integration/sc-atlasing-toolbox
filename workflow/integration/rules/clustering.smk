use rule cluster from clustering as integration_cluster with:
    input:
        zarr=rules.integration_postprocess.output.zarr
    output:
        tsv=out_dir / 'clustering' / paramspace.wildcard_pattern / 'resolutions' / '{output_type}--{resolution}.tsv',
    params:
        neighbors_key='neighbors_{output_type}',
        cluster_key_suffix='_{output_type}',
    resources:
        partition=lambda w: mcfg.get_resource(profile='gpu', resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='gpu', resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu', resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: mcfg.get_resource(profile='gpu', resource_key='gpu'),


use rule merge from clustering as integration_cluster_merge with:
    input:
        tsv=lambda wildcards: expand(
            rules.integration_cluster.output.tsv,
            resolution=mcfg.get_for_dataset(
                dataset=wildcards.dataset,
                query=['clustering', 'resolutions'],
                default=[0.4, 0.8, 1.0]
            ),
            output_type=mcfg.get_from_parameters(
                query_dict=wildcards,
                parameter_key='output_type'
            ),
            allow_missing=True
        ),
    output:
        tsv=out_dir / 'clustering' / paramspace.wildcard_pattern / 'clusters_all_resolutions.tsv',


use rule plot_umap from clustering as integration_cluster_plot_umap with:
    input:
        zarr=rules.integration_compute_umap.output.zarr,
        clusters=rules.integration_cluster_merge.output.tsv,
    output:
        png=image_dir / 'umap_clusters' / f'{paramspace.wildcard_pattern}.png',
    wildcard_constraints:
        lineage_key='((?![/]).)*',
    resources:
        partition=lambda w: mcfg.get_resource(profile='cpu', resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='cpu', resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu', resource_key='mem_mb', attempt=attempt),


rule clustering_all:
    input:
        mcfg.get_output_files(rules.integration_cluster_merge.output),
        mcfg.get_output_files(rules.integration_cluster_plot_umap.output),
