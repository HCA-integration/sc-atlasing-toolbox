module preprocessing:
    snakefile: "../../preprocessing/rules/rules.smk"
    config: config


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
        rep=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key'),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'resolutions' / '{resolution}.tsv',
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key'),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


use rule merge from clustering as clustering_merge with:
    input:
        tsv=lambda wildcards: mcfg.get_output_files(rules.clustering_cluster.output, subset_dict=dict(wildcards)),
    wildcard_constraints:
        dataset='\w+',
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'clustering.tsv'


use rule plot_umap from clustering as clustering_plot_umap with:
    input:
        zarr=rules.clustering_compute_umap.output.zarr,
        clusters=rules.clustering_merge.output.tsv,
    output:
        png=mcfg.image_dir / f'{paramspace.wildcard_pattern}.png',
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
