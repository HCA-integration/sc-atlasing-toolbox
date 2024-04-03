module plots:
    snakefile: "../../common/rules/plots.smk"
    config: config


use rule embedding from plots as label_harmonization_embedding with:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'embedding.png'
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'plot_colors']),
        basis=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'celltypist', 'use_rep']),
        ncols=2,
        wspace=0.5,
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),


use rule umap from plots as label_harmonization_umap with:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'umap.png',
        coordinates=mcfg.out_dir / paramspace.wildcard_pattern / 'celltypist' / 'label_harmonization_umap.npy',
    params:
        color=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'plot_colors']),
        use_rep=lambda w: mcfg.get_for_dataset(w.dataset, [mcfg.module_name, 'celltypist', 'use_rep']),
        ncols=2,
        wspace=0.5,
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


rule celltypist_umap:
    input:
        anndata=rules.celltypist.output.h5ad,
        group_assignment=rules.celltypist_index_reannotations.output.reannotation,
        coordinates=rules.label_harmonization_umap.output.coordinates,
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--umap.png',
        per_group=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap_per_group'),
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/celltypist_umap.py'


def get_for_organ(config, wildcards, key):
    organ = mcfg.get_from_parameters(wildcards, 'organ')
    try:
        organ_config = config['ORGANS'][organ]
    except KeyError as e:
        raise ValueError(f'Organ "{organ}" not found in config file') from e
    return organ_config[key]


rule dotplot:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
        group_assignment=rules.celltypist_index_reannotations.output.reannotation,
    output:
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--dotplot.png',
        per_group=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'dotplot_per_group'),
    params:
        marker_genes=lambda wildcards: get_for_organ(config, wildcards, 'marker_genes'),
        kwargs=dict(
            use_raw=False,
            standard_scale='var',
            # dendrogram=True,
            swap_axes=False,
        )
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'celltypist')
    script:
        '../scripts/dotplot.py'


# plot_columns = ['dataset', 'dataset_key', 'author_label_key']
# try:
#     plot_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
# except KeyError:
#     plot_datasets = []


rule dotplot_all:
    input: mcfg.get_output_files(rules.dotplot.output)
    localrule: True


rule plots_all:
    input:
        mcfg.get_output_files(rules.label_harmonization_embedding.output),
        mcfg.get_output_files(rules.label_harmonization_umap.output),
        mcfg.get_output_files(rules.celltypist_umap.output),
        mcfg.get_output_files(rules.dotplot.output),
    localrule: True
