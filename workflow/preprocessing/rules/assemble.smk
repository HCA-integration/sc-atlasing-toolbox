"""
Assemble anndata with different Preprocessing outputs
"""

use rule normalize from preprocessing as preprocessing_normalize with:
    input:
        lambda w: get_for_dataset(config, w.dataset, ['input', module_name])
    output:
        zarr=directory(out_dir / '{dataset}' / 'normalized.zarr'),
    params:
        raw_counts=lambda w: get_for_dataset(config, w.dataset, [module_name, 'raw_counts']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),


use rule highly_variable_genes from preprocessing as preprocessing_highly_variable_genes with:
    input:
        zarr=rules.preprocessing_normalize.output.zarr
    output:
        zarr=directory(out_dir / '{dataset}' / 'highly_variable_genes.zarr')
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'highly_variable_genes']),
        batch=lambda w: get_for_dataset(config, w.dataset, [module_name, 'batch']),
        lineage=lambda w: get_for_dataset(config, w.dataset, [module_name, 'lineage']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),


use rule pca from preprocessing as preprocessing_pca with:
    input:
        zarr=rules.preprocessing_highly_variable_genes.output.zarr,
        counts=rules.preprocessing_normalize.output.zarr,
    output:
        zarr=directory(out_dir / '{dataset}' / 'pca.zarr')
    params:
        scale=lambda w: get_for_dataset(config, w.dataset, [module_name, 'scale'])
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),


use rule neighbors from preprocessing as preprocessing_neighbors with:
    input:
        zarr=rules.preprocessing_pca.output.zarr
    output:
        zarr=directory(out_dir / '{dataset}' / 'neighbors.zarr')
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'neighbors'], default={}),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),


use rule umap from preprocessing as preprocessing_umap with:
    input:
        zarr=rules.preprocessing_neighbors.output.zarr,
        rep=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(out_dir / '{dataset}' / 'umap.zarr')
    params:
        neighbors_key='neighbors',
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),


def collect_files(wildcards):
    file_dict = {
        'counts': get_for_dataset(config, wildcards.dataset, ['input', module_name]),
        'normalize': rules.preprocessing_normalize.output.zarr,
        'highly_variable_genes': rules.preprocessing_highly_variable_genes.output.zarr,
        'pca': rules.preprocessing_pca.output.zarr,
        'neighbors': rules.preprocessing_neighbors.output.zarr,
        # 'umap': rules.preprocessing_umap.output.zarr,
    }
    assembly_config = get_for_dataset(config, wildcards.dataset, [module_name, 'assemble'])
    if assembly_config is None:
        return file_dict
    return {k: v for k, v in file_dict.items() if k in assembly_config}


rule assemble:
    input:
        unpack(collect_files)
    output:
        zarr=directory(out_dir / '{dataset}' / 'preprocessed.zarr')
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/assemble.py'
