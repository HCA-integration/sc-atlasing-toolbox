rule neighbors:
    input:
        h5ad=rules.run_method.output.h5ad
    output:
        h5ad=out_dir / paramspace.wildcard_pattern / 'neighbors.h5ad',
    conda:
        '../envs/scanpy_rapids.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/neighbors.py'


rule clustering:
    input:
        h5ad=rules.neighbors.output.h5ad
    output:
        tsv=out_dir / paramspace.wildcard_pattern / '_clustering' / '{resolution}.tsv',
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

        dfs = [pd.read_table(file, index_col='index') for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


rule clustering_all:
    input:clusters=expand(rules.clustering_merge.output,zip,**parameters[wildcard_names].to_dict('list'))


# per lineage

rule neighbors_per_lineage:
    input:
        h5ad=rules.run_per_lineage.output.h5ad
    output:
        h5ad=out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'neighbors.h5ad',
    conda:
        '../envs/scanpy_rapids.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/neighbors.py'


rule clustering_per_lineage:
    input:
        h5ad=rules.neighbors_per_lineage.output.h5ad
    output:
        tsv=out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / '_clustering' / '{resolution}.tsv',
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
        tsv=out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'clusters_all_resolutions.tsv',
    run:
        from functools import reduce

        dfs = [pd.read_table(file, index_col=0) for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


def collect_clustering_per_lineage(wildcards):
    return collect_lineages(
            wildcards,
            rules.clustering_per_lineage_merge.output
        )

rule clustering_per_lineage_collect:
    input: unpack(collect_clustering_per_lineage)
    output: touch(out_dir / paramspace.wildcard_pattern / 'per_lineage_clustering.done')


rule clustering_per_lineage_all:
    input:
        expand(
            rules.clustering_per_lineage_collect.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        )
