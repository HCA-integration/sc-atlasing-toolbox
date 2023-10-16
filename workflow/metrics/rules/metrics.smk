import os


# def get_integration_output(wildcards):
#     if wildcards.lineage_specific == 'global':
#         return rules.integration_run_method.output[0]
#     return rules.merge_lineage.output[0]


rule preprocess:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / paramspace_no_metric.wildcard_pattern / 'preprocessed.zarr'),
    conda:
        get_env(config, 'scanpy') # TODO use GPU accelerated neighbors
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='cpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='cpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='cpu', attempt=attempt),
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
        unintegrated=lambda wildcards: mcfg.get_for_dataset(
            dataset=wildcards.dataset,
            query=[mcfg.module_name, 'unintegrated'],
            default=rules.preprocess.output.zarr,
        ), # TODO: define optional input for metrics
        metrics_meta=workflow.source_path('../params.tsv')
    output:
        metric=mcfg.out_dir / f'{paramspace.wildcard_pattern}.tsv'
    benchmark:
        mcfg.out_dir / f'{paramspace.wildcard_pattern}.benchmark.tsv'
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch'),
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label'),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env', check_null=True),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w)),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w)),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w)),
        disk_mb=100,
        time="1-08:00:00",
    script:
        '../scripts/run.py'


rule run_all:
    input:
        mcfg.get_output_files(rules.run.output)

# rule run_per_lineage_all:
#     input:
#         expand(
#             rules.run.output,
#             zip,
#             **parameters.query('lineage_specific == "per_lineage"')[wildcard_names].to_dict('list')
#         )