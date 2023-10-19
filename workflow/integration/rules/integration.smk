use rule run_method from integration as integration_run_method with:
    message:
       """
       Integration: Run {wildcards.method} on {wildcards.dataset}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       """
    input:
        h5ad=lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'adata.zarr'),
        model=touch(directory(out_dir / paramspace.wildcard_pattern / 'model'))
    benchmark:
        out_dir / paramspace.wildcard_pattern / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts'),
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        output_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_type'),
        hyperparams=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hyperparams_dict'),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env'),
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w)),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w)),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w)),
        time="2-00:00:00",


use rule postprocess from integration as integration_postprocess with:
    input:
        zarr=rules.integration_run_method.output.zarr,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'postprocessed.zarr'),
    resources:
        partition=lambda w: mcfg.get_resource(profile='gpu', resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='gpu', resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu', resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: mcfg.get_resource(profile='gpu', resource_key='gpu'),


rule run_all:
    input:
        mcfg.get_output_files(rules.integration_postprocess.output)
