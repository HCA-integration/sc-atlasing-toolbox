rule celltypist:
    """
    Map author labels of datasets
    """
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        h5ad=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'adata.h5ad',
        reannotation=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'reannotation.tsv',
        relation=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'relation.tsv',
        model=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'model.pkl',
    params:
        author_label_key=lambda w: mcfg.get_for_dataset( w.dataset, query=[mcfg.module_name,'author_label_key']),
        dataset_key=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'dataset_key']),
        params=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'celltypist']),
        subsample=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'subsample']),
        force_scale=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'force_scale']),
    conda:
        get_env(config, 'celltypist')
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(resource_key='qos'),
        mem_mb=mcfg.get_resource(resource_key='mem_mb'),
    script:
        '../scripts/celltypist.py'


rule celltypist_index_reannotations:
    """
    Add a column for reannotations (should make relabeling easier)
    """
    input:
        reannotation=rules.celltypist.output.reannotation
    output:
        reannotation=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'reannotation_index.tsv',
    conda:
        get_env(config, 'celltypist')
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(resource_key='qos'),
        mem_mb=mcfg.get_resource(resource_key='mem_mb'),
    script:
        '../scripts/celltypist_reindex.py'


rule celltypist_plots:
    """
    Plots for celltypist output
    """
    input:
        model=rules.celltypist.output.model,
    output:
        treeplot=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--treeplot.png',
        treeplot_ordered=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--treeplot_ordered.png',
        heatmap=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--heatmap.png',
        # sankeyplot=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist' / 'sankeyplot.pdf',
    params:
        coarse_cell_type=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'author_label_key']),
    conda:
        get_env(config, 'celltypist')
    script:
        '../scripts/celltypist_plots.py'


# celltypist_columns = ['dataset', 'dataset_key', 'author_label_key', 'celltypist']
# try:
#     celltypist_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
# except KeyError:
#     celltypist_datasets = []

