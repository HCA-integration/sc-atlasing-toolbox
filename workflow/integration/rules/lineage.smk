def get_checkpoint_output(checkpoint, **kwargs):
    return Path(checkpoint.get(**kwargs).output[0])


checkpoint split_lineage:
    message:
        """
        Split lineages: Split {wildcards.dataset} by lineage key {wildcards.lineage_key}
        input: {input}
        output: {output}
        wildcards: {wildcards}
        """
    input: get_input
    output:
        directory(out_dir / 'dataset~{dataset}' / 'split_lineage,batch~{batch},lineage_key~{lineage_key}')
    params:
        label=lambda wildcards: get_params(wildcards,parameters,'label'),
        norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
    conda:
        '../envs/scanpy_rapids.yaml'
    retries: 3
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
    # shadow: 'minimal'
    script:
        '../scripts/split_anndata.py'


rule run_per_lineage:
    input:
        zarr=lambda w: get_checkpoint_output(checkpoints.split_lineage,**w) / f'{w.lineage}.zarr',
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'adata.zarr'),
        model=touch(directory(out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'model'))
    benchmark:
        out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}/benchmark.tsv'
    params:
        norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        raw_counts=lambda wildcards: get_params(wildcards,parameters,'raw_counts'),
        output_type=lambda wildcards: get_params(wildcards,parameters,'output_type'),
        hyperparams=lambda wildcards: get_params(wildcards,parameters,'hyperparams_dict'),
        env=lambda wildcards: get_params(wildcards,parameters,'env'),
    conda:
        lambda wildcards, params: f'../envs/{params.env}.yaml'
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
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'postprocessed.zarr'),
    conda:
        '../envs/scanpy_rapids.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/postprocess.py'


def collect_lineages(wildcards, pattern=rules.postprocess_per_lineage.output.zarr):
    checkpoint_output = get_checkpoint_output(checkpoints.split_lineage,**wildcards)
    lineages = glob_wildcards(str(checkpoint_output / "{lineage}.zarr")).lineage
    return {
        lineage: expand(pattern,lineage=lineage,**wildcards)
        for lineage in lineages
    }


rule merge_lineage:
    input:
        unpack(collect_lineages)
    output:
        h5mu=out_dir / paramspace.wildcard_pattern / 'lineages.h5mu',
    conda:
        '../envs/scanpy.yaml'
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
    # shadow: 'minimal'
    script:
        '../scripts/merge_anndata.py'


rule run_per_lineage_all:
    input:
        expand(
            rules.merge_lineage.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        )
