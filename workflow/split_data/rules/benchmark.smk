rule benchmark_per_dataset:
    input:
        benchmark=lambda wildcards: expand_per(rules.run_method.benchmark,parameters,wildcards,all_but(wildcard_names,'dataset')),
    output:
        benchmark=out_dir / 'dataset~{dataset}' / 'integration.benchmark.tsv'
    params:
        wildcards=lambda wildcards: parameters.query(f'dataset == "{wildcards.dataset}"')[wildcard_names]
    group:
        'integration'
    run:
        benchmark_df = pd.concat([pd.read_table(file) for file in input.benchmark]).reset_index(drop=True)
        benchmark_df = pd.concat([params.wildcards.reset_index(drop=True), benchmark_df],axis=1)
        print(benchmark_df)
        benchmark_df.to_csv(output.benchmark,sep='\t',index=False)


rule benchmark:
    input:
        benchmark=expand(rules.run_method.benchmark,zip,**parameters[wildcard_names].to_dict('list')),
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

module plots:
   snakefile: "../../common/rules/plots.smk"
   config: config


use rule barplot from plots as integration_barplot with:
    input:
        tsv=rules.benchmark.output.benchmark
    output:
        png=image_dir / 'benchmark' / '{metric}.png'
    group:
        'integration'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        #hue='hyperparams',
        facet_col='dataset',
        title='Integration runtime',
        description=wildcards_to_str,
        dodge=True,
    wildcard_constraints:
        metric='((?![/]).)*'


use rule barplot from plots as integration_barplot_per_dataset with:
    input:
        tsv=rules.benchmark_per_dataset.output.benchmark
    output:
        png=image_dir / 'benchmark' / '{dataset}' / '{metric}.png'
    group:
        'integration'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        #hue='hyperparams',
        title=lambda wildcards: f'Integration runtime for {wildcards.dataset}',
        description=wildcards_to_str,
        dodge=True,
    wildcard_constraints:
        metric='((?![/]).)*'


rule benchmark_all:
    input:
        expand(rules.integration_barplot.output,metric=['s', 'max_uss', 'mean_load']),
        expand(
            rules.integration_barplot_per_dataset.output,
            metric=['s', 'max_uss', 'mean_load'],
            dataset=parameters['dataset'].unique().tolist()
        )