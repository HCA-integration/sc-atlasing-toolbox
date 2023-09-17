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
        h5ad=lambda wildcards: module_config.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'adata.zarr'),
        model=touch(directory(out_dir / paramspace.wildcard_pattern / 'model'))
    benchmark:
        out_dir / paramspace.wildcard_pattern / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        raw_counts=lambda wildcards: get_params(wildcards,parameters,'raw_counts'),
        output_type=lambda wildcards: get_params(wildcards,parameters,'output_type'),
        hyperparams=lambda wildcards: get_params(wildcards,parameters,'hyperparams_dict'),
        env=lambda wildcards: get_params(wildcards,parameters,'env'),
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        time="2-00:00:00",


use rule postprocess from integration as integration_postprocess with:
    input:
        zarr=rules.integration_run_method.output.zarr,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'postprocessed.zarr'),
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),



rule run_all:
    input: expand(rules.integration_postprocess.output,zip,**parameters[wildcard_names].to_dict('list'))
