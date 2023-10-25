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
    script:
        '../scripts/preprocess.py'

def get_metric_input(wildcards):
    inputs = dict(
        h5mu=rules.preprocess.output.zarr,
        metrics_meta=workflow.source_path('../params.tsv')
    )
    if mcfg.get_from_parameters(wildcards, 'comparison'):
        unintegrated_file = mcfg.get_for_dataset(
            dataset=wildcards.dataset,
            query=[mcfg.module_name, 'unintegrated'],
            default=rules.preprocess.output.zarr,
        )
        inputs |= dict(unintegrated=unintegrated_file)
    return inputs


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
        unpack(get_metric_input)
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
