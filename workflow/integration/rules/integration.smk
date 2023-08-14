rule run_method:
    message:
       """
       Integration: Run {wildcards.method} on {wildcards.dataset}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       """
    input:
        h5ad=get_input
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
    conda:
        lambda wildcards, params: f'../envs/{params.env}.yaml'
    retries: 2
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        time="2-00:00:00",
    # shadow: 'minimal'
    script:
        '../scripts/methods/{wildcards.method}.py'


rule postprocess:
    input:
        zarr=rules.run_method.output.zarr,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'postprocessed.zarr'),
    conda:
        '../envs/scanpy_rapids.yaml'
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/postprocess.py'


rule run_all:
    input: expand(rules.postprocess.output,zip,**parameters[wildcard_names].to_dict('list'))
