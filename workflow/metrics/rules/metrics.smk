import os


def get_integration_output(wildcards):
    if wildcards.lineage_specific == 'global':
        return rules.integration_run_method.output[0]
    return rules.merge_lineage.output[0]


rule preprocess:
    input:
        lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'preprocessed.h5mu.zarr'),
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        time="1-00:00:00",
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
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb}
       """
    input:
        h5mu=rules.preprocess.output.zarr,
        unintegrated=lambda w: get_for_dataset(
            config,
            w.dataset,
            query=['input', 'integration'],
            default=rules.preprocess.output.zarr,
        ),
        metrics_meta=workflow.source_path('../params.tsv')
    output:
        metric=out_dir / paramspace.wildcard_pattern / '{metric}.tsv'
    params:
        env=lambda wildcards: get_params(wildcards,parameters,'env')
    conda:
        lambda wildcards, params: get_env(config, params.env)
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        mem_mb=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb'),
        disk_mb=100,
        time="1-08:00:00",
    benchmark:
        out_dir / paramspace.wildcard_pattern / '{metric}.benchmark.tsv'
    # shadow: 'copy-minimal'
    script:
        '../scripts/run.py'


rule run_all:
    input: expand(rules.run.output,zip,**parameters[wildcard_names].to_dict('list'))


# rule run_per_lineage_all:
#     input:
#         expand(
#             rules.run.output,
#             zip,
#             **parameters.query('lineage_specific == "per_lineage"')[wildcard_names].to_dict('list')
#         )