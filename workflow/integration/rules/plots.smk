######### Benchmark plots #########

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


######### UMAP and embedding plots #########

module preprocessing:
    snakefile: "../../preprocessing/rules/rules.smk"
    config: config


use rule umap from preprocessing as integration_compute_umap with:
    input:
        anndata=rules.postprocess.output.zarr,
        rep=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
    params:
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule plot_umap from preprocessing as integration_plot_umap with:
    input:
        anndata=rules.integration_compute_umap.output.zarr,
    output:
        plot=image_dir / 'umap' / f'{paramspace.wildcard_pattern}.png',
        additional_plots=directory(image_dir / 'umap' / paramspace.wildcard_pattern),
    params:
        color=lambda w: [
            get_params(w,parameters,'label'),
            get_params(w,parameters,'batch'),
            *get_for_dataset(
                config=config,
                dataset=w.dataset,
                query=[module_name, 'umap_colors'],
                default=[]
            ),
        ],
        ncols=1,
        neighbors_key=lambda w: [f'neighbors_{output_type}' for output_type in get_params(w,parameters,'output_type')],
    wildcard_constraints:
        lineage_key='((?![/]).)*',
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),


rule plots_all:
    input:
        rules.benchmark_all.input,
        expand(rules.integration_plot_umap.output,zip,**parameters[wildcard_names].to_dict('list')),


use rule umap from plots as integration_umap_lineage with:
    input:
        anndata=rules.run_per_lineage.output.zarr
    output:
        plot=image_dir / 'umap' / f'{paramspace.wildcard_pattern}' / 'lineage~{lineage}.png',
        coordinates=out_dir / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'umap_coordinates.npy'
    params:
        color=lambda wildcards: [
            get_params(wildcards,parameters,'label'),
            get_params(wildcards,parameters,'batch'),
            *get_for_dataset(
                config=config,
                dataset=wildcards.dataset,
                query = [module_name, 'umap_colors'],
                default=[]
            ),
        ],
        use_rep='X_emb',
        ncols=1,
    retries: 2
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


def collect_umap_lineages(wildcards):
    checkpoint_output = get_checkpoint_output(checkpoints.split_lineage,**wildcards)
    lineages = glob_wildcards(str(checkpoint_output / "{lineage}.zarr")).lineage
    return [
        expand(rules.integration_umap_lineage.output.plot,lineage=lineage,**wildcards)
        for lineage in lineages
    ]


rule collect_umap_lineages:
    input:
        unpack(collect_umap_lineages)
    output:
        touch(out_dir / paramspace.wildcard_pattern / 'per_lineage_umap.done')


rule plots_per_lineage_all:
    input:
        expand(
            rules.collect_umap_lineages.output,
            zip,
            **parameters.query('lineage_key != "None"')[wildcard_names].to_dict('list')
        ),