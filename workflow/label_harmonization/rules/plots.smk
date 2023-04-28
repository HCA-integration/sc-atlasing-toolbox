module plots:
   snakefile: "../../common/rules/plots.smk"
   config: config


use rule embedding from plots as label_harmonization_embedding with:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        png=out_dir / 'celltypist' / '{dataset}' / 'embedding.png'
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        basis=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule umap from plots as label_harmonization_umap with:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        png=out_dir / 'celltypist' / '{dataset}' / 'umap.png'
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        use_rep=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


plot_columns = ['dataset', 'dataset_key', 'author_label_key']
try:
    plot_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    plot_datasets = []


rule plots_all:
    input:
        expand(rules.label_harmonization_embedding.output, dataset=plot_datasets),
        expand(rules.label_harmonization_umap.output, dataset=plot_datasets)