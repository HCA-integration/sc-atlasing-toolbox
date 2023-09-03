module preprocessing:
    snakefile: "../../preprocessing/rules/rules.smk"
    config: config


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        anndata=lambda wildcards: get_for_dataset(config, wildcards.dataset, query=['input', module_name]),
        rep=lambda wildcards: get_for_dataset(config, wildcards.dataset, query=['input', module_name]),
    output:
        zarr=directory(out_dir / 'umap' / '{dataset}.zarr'),
    params:
        neighbors_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'neighbors_key']),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=lambda wildcards: get_for_dataset(config, wildcards.dataset, query=['input', module_name]),
    output:
        tsv=out_dir / 'clustering' / 'resolutions' / '{dataset}--{resolution}.tsv',
    params:
        neighbors_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'neighbors_key']),
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile='gpu',resource_key='gpu'),


rule clustering_merge:
    input:
        tsv=lambda wildcards: expand(
            rules.clustering_cluster.output.tsv,
            resolution=get_for_dataset(config, wildcards.dataset, query=[module_name, 'resolutions']),
            allow_missing=True
        ),
    output:
        tsv=out_dir / 'clustering' / '{dataset}.tsv',
    run:
        from functools import reduce

        dfs = [pd.read_table(file, index_col=0) for file in input.tsv]
        cluster_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dfs)
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


rule clustering_umap:
    input:
        zarr=rules.clustering_compute_umap.output.zarr,
        clusters=rules.clustering_merge.output.tsv,
    output:
        png=image_dir / 'umap' / '{dataset}.png',
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/plot_umap.py'
