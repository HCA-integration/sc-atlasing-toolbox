rule run:
    message:
       """
       Integration: Run {wildcards.method} on {wildcards.dataset}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu}
       """
    input:
        h5ad=get_input
    output:
        h5ad=out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/adata.h5ad',
        model=touch(directory(out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/model'))
    params:
        output_type=lambda wildcards: get_params(wildcards,parameters,'output_type'),
        hyperparams=lambda wildcards: get_params(wildcards,parameters,'hyperparams_dict'),
        env=lambda wildcards: get_params(wildcards,parameters,'env')
    conda:
        lambda wildcards, params: f'../envs/{params.env}'
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
    benchmark:
        out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/benchmark.tsv'
    shadow: 'minimal'
    script:
        '../scripts/methods/{wildcards.method}.py'


rule run_all:
    input: expand(rules.run.output,zip,**parameters[wildcard_names].to_dict('list'))
