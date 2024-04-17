class MetricNotDefinedError(RuntimeError):
    
    def __init__(self, wildcards):
        super().__init__(
            f'Metric {wildcards.metric} is not defined for dataset {wildcards.dataset} and file_id {wildcards.file_id}'
        )


rule run:
    message:
       """
       Metrics: Evaluate metric={wildcards.metric} on dataset={wildcards.dataset} and file_id={wildcards.file_id}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb}
       """
    input:
        zarr=rules.prepare.output.zarr,
        done=rules.prepare.output.done,
    output:
        metric=mcfg.out_dir / paramspace.wildcard_pattern / '{metric}.tsv'
    benchmark:
        mcfg.out_dir / paramspace.wildcard_pattern / '{metric}.benchmark.tsv'
    params:
        metric_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'metric_type', default=MetricNotDefinedError(wildcards)),
        output_types=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_types', default=MetricNotDefinedError(wildcards)),
        input_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'input_type', default=MetricNotDefinedError(wildcards)),
        comparison=lambda wildcards: mcfg.get_from_parameters(wildcards, 'comparison', default=False),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env', check_null=True, default=MetricNotDefinedError(wildcards)),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w)),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w)),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w)),
        time="1-08:00:00",
    script:
        '../scripts/run.py'


rule run_all:
    input:
        mcfg.get_output_files(rules.run.output, all_params=True)
    localrule: True
