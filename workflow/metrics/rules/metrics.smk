rule prepare:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'prepare.zarr'),
        done=touch(mcfg.out_dir / paramspace.wildcard_pattern / '.prepare.done'),
    params:
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label'),
        neighbor_args=lambda wildcards: mcfg.get_for_dataset(wildcards.dataset, ['preprocessing', 'neighbors'], default={}),
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        corrected_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'corrected', default='X'),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='gpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='gpu'),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='gpu', attempt=attempt),
        time="1-00:00:00",
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output)


# def get_metric_input(wildcards):
#     inputs = dict(
#         h5mu=rules.preprocess.output.zarr,
#         metrics_meta=workflow.source_path('../params.tsv')
#     )
#     if mcfg.get_from_parameters(wildcards, 'comparison'):
#         unintegrated_file = mcfg.get_from_parameters(query_dict=wildcards, parameter_key='unintegrated')
#         if unintegrated_file == 'None' or unintegrated_file is None:
#             wstring = ", ".join([f"{k}={v}" for k, v in wildcards.items()])
#             warnings.warn(
#                 '\nUnintegrated file is not defined for metrics module. Using default input...\n'
#                 f'wildcards: {wstring}'
#             )
#             unintegrated_file = rules.preprocess.output.zarr
#         inputs |= dict(unintegrated=unintegrated_file)
#     return inputs


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
        zarr=rules.prepare.output.zarr,
        done=rules.prepare.output.done,
    output:
        metric=mcfg.out_dir / paramspace.wildcard_pattern / '{metric}.tsv'
    benchmark:
        mcfg.out_dir / paramspace.wildcard_pattern / '{metric}.benchmark.tsv'
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch'),
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label'),
        metric_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'metric_type'),
        output_types=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_type'),
        comparison=lambda wildcards: mcfg.get_from_parameters(wildcards, 'comparison', default=False),
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
        mcfg.get_output_files(rules.run.output, all_params=True)
