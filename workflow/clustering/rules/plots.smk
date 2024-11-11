def get_cluster_keys(wildcards):
    return mcfg.get_output_files(
        '{algorithm}_{resolution}_{level}',
        subset_dict=dict(wildcards),
        all_params=True
    )


use rule plots from preprocessing as clustering_plot_umap with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}' / 'umap')
    params:
        basis='X_umap',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap_colors', default=[])
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


use rule plots from preprocessing as clustering_plot_umap_clusters with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}' / 'umap_clusters')
    params:
        basis='X_umap',
        color=get_cluster_keys,
        legend_loc='on data',
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


use rule plot_evaluation from clustering as clustering_plot_evaluation with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}' / 'evaluation')
    params:
        cluster_keys=get_cluster_keys,
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap_colors', default=[])
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: attempt * 32000,


rule plots_all:
    input:
        mcfg.get_output_files(rules.clustering_plot_umap_clusters.output),
        mcfg.get_output_files(rules.clustering_plot_umap.output),
        mcfg.get_output_files(rules.clustering_plot_evaluation.output),
    localrule: True