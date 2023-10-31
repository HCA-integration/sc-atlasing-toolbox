"""
Assemble anndata with different Preprocessing outputs
"""

use rule normalize from preprocessing as preprocessing_normalize with:
    input:
        lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'normalized.zarr'),
    params:
        raw_counts=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'raw_counts']),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),


use rule highly_variable_genes from preprocessing as preprocessing_highly_variable_genes with:
    input:
        zarr=rules.preprocessing_normalize.output.zarr
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'highly_variable_genes.zarr')
    params:
        args=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'highly_variable_genes'], default={}),
        batch=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'batch']),
        lineage=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'lineage']),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),


use rule pca from preprocessing as preprocessing_pca with:
    input:
        zarr=rules.preprocessing_highly_variable_genes.output.zarr,
        counts=rules.preprocessing_normalize.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'pca.zarr')
    params:
        args=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'pca'], default={}),
        scale=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'scale']),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),


use rule neighbors from preprocessing as preprocessing_neighbors with:
    input:
        zarr=rules.preprocessing_pca.output.zarr
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'neighbors.zarr')
    params:
        args=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'neighbors'], default={}),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),


use rule umap from preprocessing as preprocessing_umap with:
    input:
        zarr=rules.preprocessing_neighbors.output.zarr,
        rep=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'umap.zarr')
    params:
        neighbors_key='neighbors',
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),


def collect_files(wildcards):
    file_dict = {
        'counts': mcfg.get_input_file(**wildcards),
        'normalize': rules.preprocessing_normalize.output.zarr,
        'highly_variable_genes': rules.preprocessing_highly_variable_genes.output.zarr,
        'pca': rules.preprocessing_pca.output.zarr,
        'neighbors': rules.preprocessing_neighbors.output.zarr,
        'umap': rules.preprocessing_umap.output.zarr,
    }
    assembly_config = mcfg.get_for_dataset(wildcards.dataset, [mcfg.module_name, 'assemble'])
    if assembly_config is None:
        return file_dict
    return {k: v for k, v in file_dict.items() if k in assembly_config}


use rule assemble from preprocessing as preprocessing_assemble with:
    input:
        unpack(collect_files)
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'preprocessed.zarr')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=mcfg.get_resource(profile='cpu_merged',resource_key='disk_mb'),
    conda:
        get_env(config, 'scanpy')
