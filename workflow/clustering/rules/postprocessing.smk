def get_umap_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_umap', default=False):
        return rules.clustering_compute_umap.output.zarr
    return mcfg.get_input_file(**wildcards)


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        zarr=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'umap' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


use rule merge from clustering as clustering_merge with:
    input:
        zarr=get_umap_file,
        cluster_anno=lambda wildcards: mcfg.get_output_files(
            rules.clustering_cluster.output.zarr,
            subset_dict=dict(wildcards),
            all_params=True
        ),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
