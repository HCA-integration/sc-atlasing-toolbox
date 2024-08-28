def get_umap_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_umap', default=False):
        return rules.clustering_compute_umap.output.zarr
    return mcfg.get_input_file(**wildcards)


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        zarr=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'umap' / 'dataset~{dataset}' / 'file_id~{file_id}.zarr'),
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
        zarr=directory(mcfg.out_dir / 'dataset~{dataset}/file_id~{file_id}.zarr')
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),


use rule plots from preprocessing as clustering_plot_umap with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}')
    params:
        basis='X_umap',
        # color='{algorithm}_{resolution}_{level}',
        color=lambda wildcards: mcfg.get_output_files(
            '{algorithm}_{resolution}_{level}',
            subset_dict=dict(wildcards),
            all_params=True
        ),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),