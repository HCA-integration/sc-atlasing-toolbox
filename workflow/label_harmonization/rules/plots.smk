def get_neighbors_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.label_harmonization_neighbors.output.zarr
    return mcfg.get_input_file(**wildcards)


def get_umap_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_umap', default=False) \
    or mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.label_harmonization_umap.output.zarr
    return get_neighbors_file(wildcards)
    

def get_plotting_colors(wildcards):
    colors = mcfg.get_from_parameters(wildcards, 'plot_colors', default=[])
    if isinstance(colors, str):
        colors = [colors]
    return colors + ['groups', 'reannotation']


use rule neighbors from preprocessing as label_harmonization_neighbors with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'neighbors.zarr'),
    params:
        args=lambda wildcards: {
            'use_rep': mcfg.get_from_parameters(wildcards, 'cellhint').get('use_rep', 'X_pca'),
        },
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),



use rule umap from preprocessing as label_harmonization_umap with:
    input:
        anndata=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


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
        anndata=rules.label_harmonization_umap.output.zarr,
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
