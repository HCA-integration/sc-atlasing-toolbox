def get_checkpoint_output(checkpoint, **kwargs):
    return Path(checkpoint.get(**kwargs).output[0])


def get_lineages(checkpoint, **kwargs):
    checkpoint_output = get_checkpoint_output(checkpoint, **kwargs)
    return glob_wildcards(str(checkpoint_output / "lineage~{lineage}.zarr")).lineage


checkpoint split_lineage:
    message:
        """
        Split lineages: Split {wildcards.dataset} by lineage key {wildcards.lineage_key}
        input: {input}
        output: {output}
        wildcards: {wildcards}
        resources: {resources.mem_mb}MB
        """
    input: lambda wildcards: get_for_dataset(config, wildcards.dataset, query=['input', module_name])
    output:
        directory(out_dir / 'per_lineage' / 'split_lineage' / 'dataset~{dataset}--batch~{batch}--lineage_key~{lineage_key}')
    params:
        label=lambda wildcards: get_params(wildcards,parameters,'label'),
        # norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        hvg_args=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'highly_variable_genes']),
    conda:
        get_env(config, 'scanpy', gpu_env='scanpy_rapids')
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
    script:
        '../scripts/split_anndata.py'


### Preprocessing ###
pp_per_lineage_dir = out_dir / 'per_lineage' / 'preprocessing' / 'dataset~{dataset}--batch~{batch}--lineage_key~{lineage_key}={lineage}'


use rule normalize from preprocessing as preprocessing_per_normalize with:
    input:
        zarr=lambda w: get_checkpoint_output(checkpoints.split_lineage,**w) / f'lineage~{w.lineage}.zarr',
    output:
        zarr=directory(pp_per_lineage_dir / 'normalize.zarr'),
    params:
        raw_counts=lambda w: get_for_dataset(config, w.dataset, [module_name, 'raw_counts']),
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule highly_variable_genes from preprocessing as preprocessing_per_lineage_highly_variable_genes with:
    input:
        zarr=rules.preprocessing_per_normalize.output.zarr,
    output:
        zarr=directory(pp_per_lineage_dir / 'highly_variable_genes.zarr'),
    params:
        args=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'highly_variable_genes']),
        batch=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'batch']),
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule pca from preprocessing as preprocessing_per_lineage_pca with:
    input:
        zarr=rules.preprocessing_per_lineage_highly_variable_genes.output.zarr,
        counts=rules.preprocessing_per_lineage_highly_variable_genes.input.zarr,
    output:
        zarr=directory(pp_per_lineage_dir / 'pca.zarr'),
    params:
        scale=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'scale'])
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule neighbors from preprocessing as preprocessing_per_lineage_neighbors with:
    input:
        zarr=rules.preprocessing_per_lineage_pca.output.zarr
    output:
        zarr=directory(pp_per_lineage_dir / 'neighbors.zarr'),
    params:
        args=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'neighbors'], default={}),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),


use rule assemble from preprocessing as preprocessing_per_lineage_assemble with:
    input:
        unpack(
            lambda wildcards: dict(
                counts=get_checkpoint_output(checkpoints.split_lineage,**wildcards) / f'lineage~{wildcards.lineage}.zarr',
                normalize=rules.preprocessing_per_normalize.output.zarr,
                highly_variable_genes=rules.preprocessing_per_lineage_highly_variable_genes.output.zarr,
                pca=rules.preprocessing_per_lineage_pca.output.zarr,
                neighbors=rules.preprocessing_per_lineage_neighbors.output.zarr,
            )
        )
    output:
        zarr=directory(pp_per_lineage_dir / 'assembled.zarr'),
    resources:
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),


### Methods ###

rule run_per_lineage:
    input:
        zarr=rules.preprocessing_per_lineage_assemble.output.zarr,
    output:
        zarr=directory(out_dir / 'per_lineage' / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'adata.zarr'),
        model=touch(directory(out_dir / 'per_lineage' / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'model'))
    benchmark:
        out_dir / 'per_lineage' / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        raw_counts=lambda wildcards: get_params(wildcards,parameters,'raw_counts'),
        output_type=lambda wildcards: get_params(wildcards,parameters,'output_type'),
        hyperparams=lambda wildcards: get_params(wildcards,parameters,'hyperparams_dict'),
        env=lambda wildcards: get_params(wildcards,parameters,'env'),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        mem_mb=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb'),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        time="2-00:00:00",
    # shadow: 'minimal'
    script:
        '../scripts/methods/{wildcards.method}.py'


rule postprocess_per_lineage:
    input:
        zarr=rules.run_per_lineage.output.zarr
    output:
        zarr=directory(out_dir / 'per_lineage' / 'postprocess' /paramspace.wildcard_pattern / 'lineage~{lineage}' / 'postprocessed.zarr'),
    conda:
        get_env(config, 'scanpy', gpu_env='scanpy_rapids')
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/postprocess.py'


def collect_lineages(wildcards, pattern=rules.postprocess_per_lineage.output.zarr):
    return {
        lineage: expand(pattern,lineage=lineage,**wildcards)
        for lineage in get_lineages(checkpoints.split_lineage,**wildcards)
    }


rule merge_lineage:
    input:
        unpack(collect_lineages)
    output:
        zarr=directory(out_dir / 'per_lineage' / 'merge' / paramspace.wildcard_pattern / 'lineages.h5mu.zarr'),
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/merge_anndata.py'


rule run_per_lineage_all:
    input:
        expand(
            rules.merge_lineage.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        )
