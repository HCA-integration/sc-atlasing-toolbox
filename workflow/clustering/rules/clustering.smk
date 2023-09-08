module preprocessing:
    snakefile: "../../preprocessing/rules/rules.smk"
    config: config


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        anndata=lambda wildcards: get_input_file(config, wildcards, module_name),
        rep=lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        zarr=directory(out_dir / wildcard_pattern / 'umap.zarr'),
    params:
        neighbors_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'neighbors_key']),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        tsv=out_dir / wildcard_pattern / 'resolutions' / '{resolution}.tsv',
    params:
        neighbors_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'neighbors_key']),
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile='gpu',resource_key='gpu'),


use rule merge from clustering as clustering_merge with:
    input:
        tsv=lambda wildcards: expand(
            rules.clustering_cluster.output.tsv,
            resolution=get_for_dataset(config, wildcards.dataset, query=[module_name, 'resolutions']),
            allow_missing=True
        ),
    wildcard_constraints:
        dataset='\w+',
    output:
        tsv=out_dir / wildcard_pattern / 'clustering.tsv'


use rule plot_umap from clustering as clustering_plot_umap with:
    input:
        zarr=rules.clustering_compute_umap.output.zarr,
        clusters=rules.clustering_merge.output.tsv,
    output:
        png=image_dir / f'{wildcard_pattern}.png',
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
