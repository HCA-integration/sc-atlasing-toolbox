rule cellhint_plots:
    """
    CellHint plots
    * treeplot (for all and per group)
    * ordered treeplot (for all and per group)
    * heatmap
    """
    input:
        model=rules.cellhint.output.model,
    output:
        treeplot=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'treeplot'),
        treeplot_ordered=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'treeplot_ordered'),
        heatmap=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'heatmap.png',
        # sankeyplot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'sankeyplot.pdf',
    params:
        coarse_cell_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'author_label_key'),
    conda:
        get_env(config, 'cellhint')
    retries: 0
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/cellhint_plots.py'


use rule plots from preprocessing as label_harmonization_plot_umap with:
    input:
        anndata=get_umap_file,
    output:
        plots=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'umaps'),
    params:
        color=get_plotting_colors,
        basis='X_umap',
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),


rule cellhint_umap:
    input:
        anndata=get_umap_file,
        group_assignment=rules.cellhint.output.reannotation,
    output:
        plot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'umap.png',
        per_group=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'umaps_per_group'),
    params:
        # color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'plot_colors']),
        # ncols=2,
        # wspace=0.5,
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/cellhint_umap.py'


rule dotplot:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        group_assignment=rules.cellhint.output.reannotation,
    output:
        plot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'dotplot.png',
        per_group=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'dotplot_per_group'),
    params:
        marker_genes=lambda wildcards: get_marker_gene_set(mcfg, wildcards),
        kwargs=dict(
            use_raw=False,
            standard_scale='var',
            # dendrogram=True,
        )
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'cellhint')
    script:
        '../scripts/dotplot.py'


rule dotplot_all:
    input: mcfg.get_output_files(rules.dotplot.output)
    localrule: True



rule umap_all:
    input:
        mcfg.get_output_files(rules.label_harmonization_plot_umap.output),
        mcfg.get_output_files(rules.cellhint_umap.output),
    localrule: True


rule plots_all:
    input:
        mcfg.get_output_files(rules.label_harmonization_plot_umap.output),
        mcfg.get_output_files(rules.cellhint_umap.output),
        mcfg.get_output_files(rules.dotplot.output),
    localrule: True
