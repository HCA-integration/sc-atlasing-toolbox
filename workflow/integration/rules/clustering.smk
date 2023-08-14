rule clustering:
    input:
        zarr=rules.postprocess.output.zarr
    output:
        tsv=out_dir / paramspace.wildcard_pattern / 'clustering' / '{resolution}.tsv',
    conda:
        '../envs/scanpy.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/clustering.py'


rule clustering_merge:
    input:
        tsv=expand(
            rules.clustering.output.tsv,
            resolution=[0.8, 1.0, 1.4, 2.0],
            allow_missing=True
        ),
    output:
        tsv=out_dir / paramspace.wildcard_pattern / 'clusters_all_resolutions.tsv',
    run:
        from functools import reduce

        dfs = [pd.read_table(file, index_col=0) for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


rule clustering_umap:
    input:
        zarr=rules.integration_compute_umap.output.zarr,
        clusters=rules.clustering_merge.output.tsv,
    output:
        png=image_dir / 'umap_clusters' / f'{paramspace.wildcard_pattern}.png',
    conda:
        '../envs/scanpy.yaml'
    wildcard_constraints:
        lineage_key='((?![/]).)*',
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/clustering_umap.py'


rule clustering_all:
    input:
        expand(rules.clustering_merge.output,zip,**parameters[wildcard_names].to_dict('list')),
        expand(rules.clustering_umap.output,zip,**parameters[wildcard_names].to_dict('list')),


################# Per lineage clustering #################

rule clustering_per_lineage:
    input:
        zarr=rules.postprocess_per_lineage.output.zarr
    output:
        tsv=out_dir / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'clustering' / '{resolution}.tsv',
    conda:
        '../envs/scanpy.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/clustering.py'


rule clustering_per_lineage_merge:
    input:
        tsv=expand(
            rules.clustering_per_lineage.output.tsv,
            resolution=[0.8, 1.0, 1.4, 2.0],
            allow_missing=True
        ),
    output:
        tsv=out_dir / 'per_lineage' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'clusters_all_resolutions.tsv',
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
        '../envs/scanpy.yaml'
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
    output: touch(out_dir / paramspace.wildcard_pattern / 'per_lineage_clustering.done')


rule clustering_per_lineage_all:
    input:
        expand(
            rules.clustering_per_lineage_collect.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        )
