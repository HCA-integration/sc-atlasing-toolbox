rule cellhint:
    """
    Map author labels of datasets
    """
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        h5ad=mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'adata.h5ad',
        reannotation=mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'reannotation.tsv',
        relation=mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'relation.tsv',
        model=mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'model.pkl',
    params:
        author_label_key=lambda w: mcfg.get_for_dataset( w.dataset, query=[mcfg.module_name,'author_label_key']),
        dataset_key=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'dataset_key']),
        params=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'cellhint'], default={}),
        subsample=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'subsample'], default=False),
        force_scale=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'force_scale'], default=False),
    conda:
        get_env(config, 'cellhint')
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(resource_key='qos'),
        mem_mb=mcfg.get_resource(resource_key='mem_mb'),
    script:
        '../scripts/cellhint.py'


rule cellhint_index_reannotations:
    """
    Add a column for reannotations (should make relabeling easier)
    """
    input:
        reannotation=rules.cellhint.output.reannotation
    output:
        reannotation=mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'reannotation_index.tsv',
    conda:
        get_env(config, 'cellhint')
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(resource_key='qos'),
        mem_mb=mcfg.get_resource(resource_key='mem_mb'),
    script:
        '../scripts/cellhint_reindex.py'


rule cellhint_plots:
    """
    Plots for cellhint output
    """
    input:
        model=rules.cellhint.output.model,
    output:
        treeplot=mcfg.image_dir / paramspace.wildcard_pattern / 'cellhint--treeplot.png',
        treeplot_ordered=mcfg.image_dir / paramspace.wildcard_pattern / 'cellhint--treeplot_ordered.png',
        heatmap=mcfg.image_dir / paramspace.wildcard_pattern / 'cellhint--heatmap.png',
        # sankeyplot=mcfg.image_dir / paramspace.wildcard_pattern / 'cellhint' / 'sankeyplot.pdf',
    params:
        coarse_cell_type=lambda w: mcfg.get_for_dataset(w.dataset, query=[mcfg.module_name,'author_label_key']),
    conda:
        get_env(config, 'cellhint')
    script:
        '../scripts/cellhint_plots.py'


# cellhint_columns = ['dataset', 'dataset_key', 'author_label_key', 'cellhint']
# try:
#     cellhint_datasets = parameters[cellhint_columns].dropna()['dataset'].unique()
# except KeyError:
#     cellhint_datasets = []

