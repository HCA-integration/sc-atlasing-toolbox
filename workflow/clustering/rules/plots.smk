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
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]),
        # outlier_factor=10,
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
        # outlier_factor=10,
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


def get_mem_mb(attempt, profile):
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    return min(100_000, int(mem_mb // 5))


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
    threads:
        lambda wildcards: max(1, min(10, len(mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]))))
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt=attempt, profile='cpu'),

rule plots_all:
    input:
        mcfg.get_output_files(rules.clustering_plot_umap_clusters.output),
        mcfg.get_output_files(rules.clustering_plot_umap.output),
        mcfg.get_output_files(rules.clustering_plot_evaluation.output),
    localrule: True