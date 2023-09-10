rule celltypist:
    """
    Map author labels of datasets
    """
    input:
        anndata=lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        h5ad=out_dir / wildcard_pattern / 'celltypist' / 'adata.h5ad',
        reannotation=out_dir / wildcard_pattern / 'celltypist' / 'reannotation.tsv',
        relation=out_dir / wildcard_pattern / 'celltypist' / 'relation.tsv',
        model=out_dir / wildcard_pattern / 'celltypist' / 'model.pkl',
    params:
        author_label_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'author_label_key']),
        dataset_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'dataset_key']),
        params=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'celltypist']),
        subsample=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'subsample']),
        force_scale=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'force_scale']),
    conda:
        get_env(config, 'celltypist')
    resources:
        partition=lambda w: get_resource(config,resource_key='partition'),
        qos=lambda w: get_resource(config,resource_key='qos'),
        mem_mb=lambda w: get_resource(config,resource_key='mem_mb'),
    script:
        '../scripts/celltypist.py'


rule celltypist_index_reannotations:
    """
    Add a column for reannotations (should make relabeling easier)
    """
    input:
        reannotation=rules.celltypist.output.reannotation
    output:
        reannotation=out_dir / wildcard_pattern / 'celltypist' / 'reannotation_index.tsv',
    conda:
        get_env(config, 'celltypist')
    resources:
        partition=lambda w: get_resource(config,resource_key='partition'),
        qos=lambda w: get_resource(config,resource_key='qos'),
        mem_mb=lambda w: get_resource(config,resource_key='mem_mb'),
    script:
        '../scripts/celltypist_reindex.py'


rule celltypist_plots:
    """
    Plots for celltypist output
    """
    input:
        model=rules.celltypist.output.model,
    output:
        treeplot=image_dir / wildcard_pattern / 'celltypist--treeplot.png',
        treeplot_ordered=image_dir / wildcard_pattern / 'celltypist--treeplot_ordered.png',
        heatmap=image_dir / wildcard_pattern / 'celltypist--heatmap.png',
        # sankeyplot=out_dir / wildcard_pattern / 'celltypist' / 'sankeyplot.pdf',
    params:
        coarse_cell_type=lambda w: get_for_dataset(config, w.dataset, query=[module_name,'author_label_key']),
    conda:
        get_env(config, 'celltypist')
    script:
        '../scripts/celltypist_plots.py'


celltypist_columns = ['dataset', 'dataset_key', 'author_label_key', 'celltypist']
try:
    celltypist_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    celltypist_datasets = []

