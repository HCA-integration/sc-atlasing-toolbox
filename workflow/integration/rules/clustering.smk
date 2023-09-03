use rule cluster from clustering as integration_cluster with:
    input:
        zarr=rules.integration_postprocess.output.zarr
    output:
        tsv=out_dir / 'clustering' / paramspace.wildcard_pattern / 'resolutions' / '{output_type}--{resolution}.tsv',
    params:
        neighbors_key='neighbors_{output_type}',
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile='gpu',resource_key='gpu'),


use rule merge from clustering as integration_cluster_merge with:
    input:
        tsv=lambda wildcards: expand(
            rules.integration_cluster.output.tsv,
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
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


rule clustering_all:
    input:
        expand(rules.integration_cluster_merge.output,zip,**parameters[wildcard_names].to_dict('list')),
        expand(rules.integration_cluster_plot_umap.output,zip,**parameters[wildcard_names].to_dict('list')),


################# Per lineage clustering #################

rule clustering_per_lineage:
    input:
        zarr=rules.postprocess_per_lineage.output.zarr
    output:
        tsv=out_dir / 'clustering' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'resolutions' / '{resolution}.tsv',
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile='gpu',resource_key='gpu'),
    script:
        '../scripts/clustering.py'


rule clustering_per_lineage_merge:
    input:
        tsv=expand(
            rules.clustering_per_lineage.output.tsv,
            resolution=[0.4, 0.8, 1.0],
            allow_missing=True
        ),
    output:
        tsv=out_dir / 'clustering' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'clusters_all_resolutions.tsv',
    run:
        from functools import reduce

        dfs = [pd.read_table(file, index_col=0) for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')

rule clustering_per_lineage_umap:
    input:
        zarr=rules.integration_compute_umap_lineage.output.zarr,
        clusters=rules.clustering_per_lineage_merge.output.tsv,
    output:
        png=image_dir / 'umap_clusters' / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}.png',
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/clustering_umap.py'


rule clustering_per_lineage_collect:
    input:
        unpack(lambda w: collect_lineages(w, rules.clustering_per_lineage_merge.output)),
        unpack(lambda w: collect_lineages(w, rules.clustering_per_lineage_umap.output)),
    output: touch(image_dir / 'umap_clusters' / 'per_lineage' / paramspace.wildcard_pattern / 'clustering.done')


rule clustering_per_lineage_all:
    input:
        expand(
            rules.clustering_per_lineage_collect.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        )
