module preprocessing:
    snakefile: "../../preprocessing/rules/rules.smk"
    config: config


def get_neighbors_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.clustering_compute_neighbors.output.zarr
    return mcfg.get_input_file(**wildcards)


def get_umap_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'umap_exists', default=False):
        return mcfg.get_input_file(**wildcards)
    return rules.clustering_compute_umap.output.zarr


use rule neighbors from preprocessing as clustering_compute_neighbors with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'neighbors' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=get_neighbors_file,
    output:
        parquet=mcfg.out_dir / 'resolutions' / paramspace.wildcard_pattern / '{resolution}.parquet',
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
        algorithm=lambda wildcards: mcfg.get_from_parameters(wildcards, 'algorithm', default='leiden'),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),



use rule umap from preprocessing as clustering_compute_umap with:
    input:
        zarr=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'umap' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),


use rule merge from clustering as clustering_merge with:
    input:
        zarr=get_umap_file, # rules.clustering_compute_umap.output.zarr,
        parquet=lambda wildcards: mcfg.get_output_files(
            rules.clustering_cluster.output.parquet,
            subset_dict=dict(wildcards),
            all_params=True,
        ),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr')
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
        algorithm=lambda wildcards: mcfg.get_from_parameters(wildcards, 'algorithm', default='leiden')


use rule plots from preprocessing as clustering_plot_umap with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap')
    params:
        basis='X_umap',
        color=lambda wildcards: [
            f"{mcfg.get_from_parameters(wildcards, 'algorithm', default='leiden')}_{resolution}"
            for resolution in mcfg.get_from_parameters(wildcards, 'resolution', default=[], single_value=False)
        ],
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
