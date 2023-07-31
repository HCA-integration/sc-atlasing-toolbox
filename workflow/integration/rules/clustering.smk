rule clustering:
    input:
        h5ad=rules.run_method.output.h5ad
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
        tsv=expand(rules.clustering.output.tsv, resolution=[0.8, 1.0, 1.4, 2.0], allow_missing=True),
    output:
        tsv=out_dir / paramspace.wildcard_pattern / 'clusters_all_resolutions.tsv',
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    run:
        from functools import reduce

        dfs = [pd.read_table(file, index_col='index') for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


rule clustering_all:
    input:clusters=expand(rules.clustering_merge.output,zip,**parameters[wildcard_names].to_dict('list'))