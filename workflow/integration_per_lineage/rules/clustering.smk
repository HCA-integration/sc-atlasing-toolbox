use rule cluster from clustering as integration_per_lineage_cluster with:
    input:
        zarr=rules.integration_per_lineage_postprocess.output.zarr
    output:
        tsv=out_dir / 'clustering' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'resolutions' / '{output_type}--{resolution}.tsv',
    params:
        neighbors_key='neighbors_{output_type}',
        cluster_key_suffix='_{output_type}',
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile='gpu',resource_key='gpu'),


use rule merge from clustering as integration_per_lineage_cluster_merge with:
    input:
        tsv=lambda wildcards: expand(
            rules.integration_per_lineage_cluster.output.tsv,
            resolution=get_for_dataset(
                config,
                wildcards.dataset,
                query=['clustering', 'resolutions'],
                default=[0.4, 0.8, 1.0]
            ),
            output_type=get_params(wildcards, parameters, 'output_type'),
            allow_missing=True
        ),
    output:
        tsv=out_dir / 'clustering' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'clusters_all_resolutions.tsv',


use rule plot_umap from clustering as integration_per_lineage_cluster_plot_umap with:
    input:
        zarr=rules.integration_per_lineage_compute_umap.output.zarr,
        clusters=rules.integration_per_lineage_cluster_merge.output.tsv,
    output:
        png=image_dir / 'umap_clusters' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}.png',
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


rule cluster_collect:
    input:
        unpack(lambda w: collect_lineages(w, rules.integration_per_lineage_cluster_merge.output)),
        unpack(lambda w: collect_lineages(w, rules.integration_per_lineage_cluster_plot_umap.output)),
    output: touch(image_dir / 'umap_clusters' / 'per_lineage' / paramspace.wildcard_pattern / 'clustering.done')


rule cluster_all:
    input:
        expand(rules.cluster_collect.output, zip, **parameters[wildcard_names].to_dict('list'))
