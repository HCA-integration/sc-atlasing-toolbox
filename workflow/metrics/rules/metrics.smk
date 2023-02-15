import os


def get_integration_output(wildcards):
    if wildcards.lineage_specific == 'true':
        return rules.merge_lineage.output.h5ad
    return rules.integration_run.output.h5ad


# TODO: lineage-specific?
rule preprocess:
    input:
        h5ad=get_integration_output
    output:
        h5mu=out_dir / paramspace.wildcard_pattern / 'preprocessed.h5mu'
    conda:
        lambda wildcards, params: f'../envs/scib.yaml'
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
    # shadow: 'copy-minimal'
    script:
        '../scripts/preprocess.py'


rule run:
    message:
       """
       Metrics: Evaluate {wildcards.metric} on {wildcards.dataset}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu}
       """
    input:
        h5mu=rules.preprocess.output.h5mu,
        metrics_meta=workflow.source_path('../params.tsv')
    output:
        metric=out_dir / paramspace.wildcard_pattern / '{metric}.tsv'
    params:
        env=lambda wildcards: get_params(wildcards,parameters,'env')
    conda:
        lambda wildcards, params: f'../envs/{params.env}.yaml'
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        mem_mb=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb'),
        disk_mb=100
    benchmark:
        out_dir / paramspace.wildcard_pattern / '{metric}.benchmark.tsv'
    # shadow: 'copy-minimal'
    script:
        '../scripts/run.py'
