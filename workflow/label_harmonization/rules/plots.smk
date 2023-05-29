module plots:
   snakefile: "../../common/rules/plots.smk"
   config: config


use rule embedding from plots as label_harmonization_embedding with:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        plot=image_dir / '{dataset}' / 'embedding.png'
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        basis=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
        ncols=2,
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule umap from plots as label_harmonization_umap with:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        plot=image_dir / '{dataset}' / 'umap.png',
        coordinates=out_dir / 'celltypist' / '{dataset}' / 'label_harmonization_umap.npy',
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        use_rep=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
        ncols=2,
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


use rule umap from plots as label_harmonization_celltypist_umap with:
    input:
        anndata=rules.celltypist.output.h5ad,
    output:
        plot=image_dir / '{dataset}' / 'celltypist--umap.png',
        coordinates=out_dir / 'celltypist' / '{dataset}' / 'label_harmonization_celltypist_umap.npy',
    params:
        color='high_hierarchy',
        use_rep=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
        ncols=1,
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


rule dotplot:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
        group_assignment=rules.celltypist.output.reannotation,
    output:
        plot=image_dir / '{dataset}' / 'celltypist--dotplot.png',
    params:
        marker_genes=lambda w: config['ORGANS'][get_for_dataset(config, w.dataset, [module_name, 'organ'])]['marker_genes'],
        kwargs=dict(
            use_raw=False,
            standard_scale='var',
            dendrogram=False,
            swap_axes=True,
        )
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/dotplot.py'


plot_columns = ['dataset', 'dataset_key', 'author_label_key']
try:
    plot_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    plot_datasets = []


rule plots_all:
    input:
        expand(rules.label_harmonization_embedding.output, dataset=plot_datasets),
        expand(rules.label_harmonization_umap.output, dataset=plot_datasets),
        expand(rules.label_harmonization_celltypist_umap.output, dataset=plot_datasets),
        expand(rules.dotplot.output, dataset=plot_datasets),
