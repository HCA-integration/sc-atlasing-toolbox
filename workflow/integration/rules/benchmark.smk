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


rule plots:
    input: expand(rules.barplot.output,metric=['s', 'max_uss', 'mean_load'])
