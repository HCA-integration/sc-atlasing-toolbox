"""
Assemble anndata with different Preprocessing outputs
"""

use rule normalize from preprocessing as preprocessing_normalize with:
    input:
        lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'normalized.zarr'),
    params:
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        backed=lambda wildcards: mcfg.get_from_parameters(wildcards, 'backed'),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask'),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt, factor=1),


use rule highly_variable_genes from preprocessing as preprocessing_highly_variable_genes with:
    input:
        zarr=rules.preprocessing_normalize.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'highly_variable_genes.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'highly_variable_genes', default={}),
        backed=lambda wildcards: mcfg.get_from_parameters(wildcards, 'backed'),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask'),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt, factor=1),


use rule extra_hvgs from preprocessing as preprocessing_extra_hvgs with:
    input:
        zarr=rules.preprocessing_normalize.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'extra_hvgs.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'highly_variable_genes', default={}),
        extra_hvgs=lambda wildcards: mcfg.get_from_parameters(wildcards, 'extra_hvgs', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt, factor=1),


use rule pca from preprocessing as preprocessing_pca with:
    input:
        zarr=rules.preprocessing_highly_variable_genes.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'pca.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca', default={}),
        scale=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scale', default=False),
        backed=lambda wildcards: mcfg.get_from_parameters(wildcards, 'backed'),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask'),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt, factor=1),


use rule neighbors from preprocessing as preprocessing_neighbors with:
    input:
        zarr=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'neighbors.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),


use rule umap from preprocessing as preprocessing_umap with:
    input:
        zarr=rules.preprocessing_neighbors.output.zarr,
        rep=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
        done=touch(mcfg.out_dir / paramspace.wildcard_pattern / 'umap.done'),
    params:
        # args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap', default={}),  # TODO use args instead of direct params
        neighbors_key='neighbors',
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt),


def collect_files(wildcards):
    file_dict = {
        # 'counts': mcfg.get_input_file(**wildcards),
        'normalize': rules.preprocessing_normalize.output.zarr,
        'highly_variable_genes': rules.preprocessing_highly_variable_genes.output.zarr,
        'extra_hvgs': rules.preprocessing_extra_hvgs.output.zarr,
        'pca': rules.preprocessing_pca.output.zarr,
        'neighbors': rules.preprocessing_neighbors.output.zarr,
        'umap': rules.preprocessing_umap.output.zarr,
    }
    assembly_config = mcfg.get_from_parameters(wildcards, 'assemble')
    if assembly_config is None:
        return file_dict
    return {
        k: v for k, v in file_dict.items()
        if k in assembly_config
    }


use rule assemble from preprocessing as preprocessing_assemble with:
    input:
        unpack(collect_files)
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'preprocessed.zarr')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    retries: 1
    conda:
        get_env(config, 'scanpy')
