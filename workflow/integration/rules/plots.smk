module plots:
   snakefile: "../../common/rules/plots.smk"
   config: config


rule benchmark:
    input:
        benchmark=expand(rules.run.benchmark,zip,**parameters[wildcard_names].to_dict('list'))
    output:
        benchmark=out_dir / 'integration.benchmark.tsv'
    params:
        wildcards=parameters[wildcard_names]
    group:
        'integration'
    run:
        benchmark_df = pd.concat([pd.read_table(file) for file in input.benchmark]).reset_index(drop=True)
        benchmark_df = pd.concat([params.wildcards.reset_index(drop=True), benchmark_df],axis=1)
        print(benchmark_df)
        benchmark_df.to_csv(output.benchmark,sep='\t',index=False)


use rule barplot from plots as integration_barplot with:
    input:
        tsv=rules.benchmark.output.benchmark
    output:
        png=out_dir / 'plots' / '{metric}.png'
    group:
        'integration'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        #hue='hyperparams',
        facet_col='dataset',
        title='Integration benchmark methods',
        description=wildcards_to_str,
        dodge=True,


use rule umap from plots as integration_umap with:
    input:
        anndata=rules.run.output.h5ad
    output:
        png=out_dir / paramspace.wildcard_pattern / '{method}_umap.png'
    params:
        color=lambda wildcards: [
            get_params(wildcards,parameters,'label'),
            get_params(wildcards,parameters,'batch')
        ],
        use_rep='X_emb'
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


rule plots_all:
    input:
        expand(rules.integration_barplot.output,metric=['s', 'max_uss', 'mean_load']),
        expand(rules.integration_umap.output,zip,**parameters[wildcard_names].to_dict('list')),


use rule umap from plots as integration_umap_lineage with:
    input:
        anndata=rules.run_per_lineage.output.h5ad
    output:
        png=out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / '{method}-{lineage}_umap.png'
    params:
        color=lambda wildcards: [
            get_params(wildcards,parameters,'label'),
            get_params(wildcards,parameters,'batch')
        ],
        use_rep='X_emb'
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


def collect_umap_lineages(wildcards):
    checkpoint_output = get_checkpoint_output(checkpoints.split_lineage,**wildcards)
    lineages = glob_wildcards(str(checkpoint_output / "{lineage}.h5ad")).lineage
    return [
        expand(rules.integration_umap_lineage.output.png,lineage=lineage,**wildcards)
        for lineage in lineages
    ]


rule collect_umap_lineages:
    input:
        unpack(collect_umap_lineages)
    output:
        touch(out_dir / paramspace.wildcard_pattern / 'umap_lineages.done')


rule plots_per_lineage_all:
    input:
        expand(
            rules.collect_umap_lineages.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        ),