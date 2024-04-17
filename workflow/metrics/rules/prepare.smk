rule prepare:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'prepare.zarr'),
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
    localrule: True


use rule cluster from clustering as metrics_cluster with:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        tsv=mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'cluster_resolutions' / '{resolution}.tsv',
    params:
        algorithm=lambda wildcards: mcfg.get_from_parameters(
            wildcards,
            'cluster_algorithm',
            default='leiden',
            wildcards_sub=list(wildcards.keys())
        ),
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='gpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='gpu'),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='gpu', attempt=attempt),


use rule merge from clustering as metrics_cluster_collect with:
    input:
        zarr=rules.prepare.output.zarr,
        tsv=expand(
            rules.metrics_cluster.output.tsv,
            resolution=[2 * (x + 1) / 10 for x in range(10)],
            allow_missing=True,
        ),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'clustered.zarr'),
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='gpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='gpu'),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='gpu', attempt=attempt),


rule cluster_all:
    input:
        mcfg.get_output_files(rules.metrics_cluster_collect.output)
    localrule: True