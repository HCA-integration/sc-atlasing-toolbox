module plots:
    snakefile: "../../common/rules/plots.smk"
    config: config


use rule embedding from plots as label_harmonization_embedding with:
    input:
        anndata=lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        plot=image_dir / wildcard_pattern / 'embedding.png'
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        basis=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
        ncols=2,
        wspace=0.5,
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


use rule umap from plots as label_harmonization_umap with:
    input:
        anndata=lambda wildcards: get_input_file(config, wildcards, module_name),
    output:
        plot=image_dir / wildcard_pattern / 'umap.png',
        coordinates=out_dir / wildcard_pattern / 'celltypist' / 'label_harmonization_umap.npy',
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'plot_colors']),
        use_rep=lambda w: get_for_dataset(config, w.dataset, [module_name, 'celltypist', 'use_rep']),
        ncols=2,
        wspace=0.5,
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),


rule celltypist_umap:
    input:
        anndata=rules.celltypist.output.h5ad,
        group_assignment=rules.celltypist_index_reannotations.output.reannotation,
        coordinates=rules.label_harmonization_umap.output.coordinates,
    output:
        plot=image_dir / wildcard_pattern / 'celltypist--umap.png',
        per_group=directory(image_dir / wildcard_pattern / 'umap_per_group'),
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/celltypist_umap.py'


def get_for_organ(config, wildcards, module_name, key):
    organ = get_for_dataset(config, wildcards.dataset, [module_name, 'organ'])
    try:
        organ_config = config['ORGANS'][organ]
    except KeyError as e:
        raise ValueError(f'Organ "{organ}" not found in config file') from e
    return organ_config[key]


rule dotplot:
    input:
        anndata=lambda wildcards: get_input_file(config, wildcards, module_name),
        group_assignment=rules.celltypist_index_reannotations.output.reannotation,
    output:
        plot=image_dir / wildcard_pattern / 'celltypist--dotplot.png',
        per_group=directory(image_dir / wildcard_pattern / 'dotplot_per_group'),
    params:
        marker_genes=lambda wildcards: get_for_organ(config, wildcards, module_name, 'marker_genes'),
        kwargs=dict(
            use_raw=False,
            standard_scale='var',
            # dendrogram=True,
            swap_axes=False,
        )
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
    conda:
        get_env(config, 'celltypist')
    script:
        '../scripts/dotplot.py'


plot_columns = ['dataset', 'dataset_key', 'author_label_key']
try:
    plot_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    plot_datasets = []


rule dotplot_all:
    input: expand(rules.dotplot.output, zip, **input_files)


rule plots_all:
    input:
        expand(rules.label_harmonization_embedding.output, zip, **input_files),
        expand(rules.label_harmonization_umap.output, zip, **input_files),
        expand(rules.celltypist_umap.output, zip, **input_files),
        expand(rules.dotplot.output, zip, **input_files),
