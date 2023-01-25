rule preprocess:
    input:
        h5ad=rules.integration_run.output.h5ad,
    output:
        h5ad=out_dir / 'datasets/{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/preprocessed.h5ad'
    conda:
        lambda wildcards, params: f'../envs/scanpy.yaml'
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
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
        h5ad=rules.preprocess.output.h5ad,
        metrics_meta=workflow.source_path('../params.tsv')
    output:
        metric=out_dir / 'datasets/{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/{metric}.tsv'
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
        out_dir / 'datasets/{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/{metric}.benchmark.tsv'
    script:
        '../scripts/metrics/{wildcards.metric}.py'
