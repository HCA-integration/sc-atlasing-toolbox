def get_plotting_colors(wildcards):
    colors = mcfg.get_from_parameters(wildcards, 'plot_colors', default=[])
    if isinstance(colors, str):
        colors = [colors]
    return colors + ['groups', 'reannotation']


use rule plots from preprocessing as label_harmonization_plot_umap with:
    input:
        anndata=get_umap_file,
    output:
        plots=directory(image_dir / paramspace.wildcard_pattern / 'umaps'),
    params:
        color=get_plotting_colors,
        basis='X_umap',
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),



checkpoint split_cellhint_groups:
    input:
        zarr=rules.cellhint.output.zarr,
    output:
        group_setup=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'groups'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/split_cellhint_groups.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.split_cellhint_groups.get(**wildcards).output[0]}/{{group}}.yaml'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        group=glob_wildcards(checkpoint_output).group,
        allow_missing=True
    )


rule cellhint_umap_per_group:
    input:
        anndata=get_umap_file,
        group_assignment=rules.cellhint.output.reannotation,
        group=get_checkpoint_output,
        splits=checkpoints.split_cellhint_groups.rule.output,
    output:
        png=directory(image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'umaps'),
    params:
        author_label_key=lambda w: mcfg.get_from_parameters({k: w[k] for k in ('dataset', 'file_id')}, 'author_label_key'),
        dataset_key=lambda w: mcfg.get_from_parameters({k: w[k] for k in ('dataset', 'file_id')}, 'dataset_key'),
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/cellhint_umap.py'


rule cellhint_plots:
    """
    CellHint plots
    * treeplot (for all and per group)
    * ordered treeplot (for all and per group)
    * heatmap
    """
    input:
        model=rules.cellhint.output.model,
        group=get_checkpoint_output,
        splits=checkpoints.split_cellhint_groups.rule.output,
    output:
        treeplot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'treeplot.png',
        treeplot_ordered=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'treeplot_ordered.png',
        heatmap=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'heatmap.png',
        # sankeyplot=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'sankeyplot.pdf',
    conda:
        get_env(config, 'cellhint')
    retries: 0
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/cellhint_plots.py'



rule cellhint_dotplot:
    input:
        zarr=rules.cellhint.output.zarr,
        group=get_checkpoint_output,
        splits=checkpoints.split_cellhint_groups.rule.output,
    output:
        png=image_dir / paramspace.wildcard_pattern / 'cellhint' / 'group~{group}' / 'dotplot.png',
    params:
        author_label_key=lambda w: mcfg.get_from_parameters({k: w[k] for k in ('dataset', 'file_id')}, 'author_label_key'),
        marker_genes=lambda wildcards: get_marker_gene_set(
            mcfg,
            wildcards={k: wildcards[k] for k in ('dataset', 'file_id')},
            flatten=True,
        ),
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


rule collect_plots:
    input:
        lambda wildcards: get_from_checkpoint(wildcards, rules.cellhint_plots.output),
        lambda wildcards: get_from_checkpoint(wildcards, rules.cellhint_umap_per_group.output),
        lambda wildcards: get_from_checkpoint(wildcards, rules.cellhint_dotplot.output),
    output:
        touch(mcfg.out_dir / paramspace.wildcard_pattern / 'cellhint' / 'plots.done'),
    localrule: True


rule plots_all:
    input:
        mcfg.get_output_files(rules.collect_plots.output),
        mcfg.get_output_files(rules.label_harmonization_plot_umap.output),
    localrule: True
