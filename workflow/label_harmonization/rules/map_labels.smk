rule celltypist:
    """
    Map author labels of datasets
    """
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name]),
    output:
        h5ad=out_dir / 'celltypist' / '{dataset}' / 'adata.h5ad',
        reannotation=out_dir / 'celltypist' / '{dataset}' / 'reannotation.tsv',
        relation=out_dir / 'celltypist' / '{dataset}' / 'relation.tsv',
        model=out_dir / 'celltypist' / '{dataset}' / 'model.pkl',
    params:
        author_label_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'author_label_key']),
        dataset_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'dataset_key']),
        params=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'celltypist']),
        subsample=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'subsample']),
        force_scale=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'force_scale']),
    conda:
        '../envs/celltypist.yaml'
    resources:
        partition=lambda w: get_resource(config,resource_key='partition'),
        qos=lambda w: get_resource(config,resource_key='qos'),
        mem_mb=lambda w: get_resource(config,resource_key='mem_mb'),
    script:
        '../scripts/celltypist.py'


rule celltypist_plots:
    """
    Plots for celltypist output
    """
    input:
        model=rules.celltypist.output.model,
    output:
        treeplot=image_dir / '{dataset}' / 'celltypist--treeplot.png',
        treeplot_ordered=image_dir / '{dataset}' / 'celltypist--treeplot_ordered.png',
        heatmap=image_dir / '{dataset}' / 'celltypist--heatmap.png',
        # sankeyplot=out_dir / 'celltypist' / '{dataset}_sankeyplot.pdf',
    params:
        coarse_cell_type=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'author_label_key']),
    conda:
        '../envs/celltypist.yaml'
    script:
        '../scripts/celltypist_plots.py'


celltypist_columns = ['dataset', 'dataset_key', 'author_label_key', 'celltypist']
try:
    celltypist_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    celltypist_datasets = []

rule celltypist_all:
    input:
        expand(rules.celltypist_plots.output, dataset=celltypist_datasets)
