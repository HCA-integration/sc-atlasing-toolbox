def get_checkpoint_output(checkpoint, **kwargs):
    return Path(checkpoint.get(**kwargs).output[0])


checkpoint split_lineage:
    message:
        """
        {wildcards}
        {params}
        """
    input:
        h5ad=get_input
    output:
        directory(out_dir / '{dataset}' / 'lineages')
    params:
        split_key=lambda wildcards: get_params(wildcards,parameters,'lineage_key')
    conda:
        '../envs/scanpy.yaml'
    shadow: 'minimal'
    script:
        '../scripts/split_anndata.py'


rule run_per_lineage:
    input:
        h5ad=lambda w: get_checkpoint_output(checkpoints.split_lineage,dataset=w.dataset) / f'{w.lineage}.h5ad',
    output:
        h5ad=out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/lineages/{lineage}/adata.h5ad',
        model=touch(
            directory(
                out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/lineages/{lineage}/model'
            )
        )
    params:
        output_type=lambda w: get_params(w,parameters,column='output_type'),
        hyperparams=lambda w: get_params(w,parameters,column='hyperparams_dict'),
        env=lambda w: get_params(w,parameters,column='env')
    conda:
        lambda wildcards, params: f'../envs/{params.env}'
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
    benchmark:
        out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/lineages/{lineage}/benchmark.tsv'
    shadow: 'minimal'
    script:
        '../scripts/methods/{wildcards.method}.py'


def collect_lineages(wildcards):
    checkpoint_output = get_checkpoint_output(checkpoints.split_lineage,**wildcards)
    targets = expand(
        rules.run_per_lineage.output.h5ad,
        lineage=glob_wildcards(str(checkpoint_output / "{lineage}.h5ad")).lineage,
        **wildcards,
    )
    return targets


rule merge_lineage:
    input:
        collect_lineages
    output:
        h5ad=out_dir / '{dataset}/{method}/batch={batch},label={label},hyperparams={hyperparams}/lineages.h5ad',
    conda:
        '../envs/scanpy.yaml'
    shadow: 'minimal'
    script:
        '../scripts/merge_anndata.py'


rule run_lineages_all:
    input: expand(rules.merge_lineage.output,zip,**parameters[wildcard_names].to_dict('list'))
