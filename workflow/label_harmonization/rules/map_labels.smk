rule cellhint:
    """
    Map author labels of datasets
    """
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'cellhint' / 'adata.zarr'),
        reannotation=out_dir / paramspace.wildcard_pattern / 'cellhint' / 'reannotation.tsv',
        relation=out_dir / paramspace.wildcard_pattern / 'cellhint' / 'relation.tsv',
        model=out_dir / paramspace.wildcard_pattern / 'cellhint' / 'model.pkl',
    params:
        author_label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'author_label_key'),
        dataset_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dataset_key'),
        params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cellhint', default={}),
        subsample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'subsample', default=False),
        force_scale=lambda wildcards: mcfg.get_from_parameters(wildcards, 'force_scale', default=False),
    conda:
        get_env(config, 'cellhint')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb',attempt=attempt),
    script:
        '../scripts/cellhint.py'


rule cellhint_plots:
    """
    Plots for cellhint output
    """
    input:
        model=rules.cellhint.output.model,
    output:
        treeplot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'treeplot.png',
        treeplot_ordered=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'treeplot_ordered.png',
        heatmap=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'heatmap.png',
        # sankeyplot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'sankeyplot.pdf',
    params:
        coarse_cell_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'author_label_key'),
    conda:
        get_env(config, 'cellhint')
    script:
        '../scripts/cellhint_plots.py'
