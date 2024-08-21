rule rank_genes_groups:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'groups' / paramspace.wildcard_pattern / 'group={group}.zarr'),
        tsv=mcfg.out_dir / 'groups' / paramspace.wildcard_pattern / 'group={group}--marker_genes.tsv',
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'rank_genes_groups', default={}),
        sample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
    conda:
        get_env(config, 'scanpy') #, gpu_env='rapids_singlecell')
    resources:
        # partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        # qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        # gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/rank_genes_groups.py'


rule plot_user:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        dotplot=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'user_markers'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot', default={}),
        markers=lambda wildcards: get_marker_gene_set(mcfg, wildcards),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/plot_user.py'


rule plot:
    input:
        zarr=rules.rank_genes_groups.output.zarr,
    output:
        rankplot=mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'rank_plot.png',
        dotplot=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'dotplot'),
        matrixplot=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'matrixplot'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot', default={}),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/plot.py'


rule collect:
    input:
        lambda wildcards: mcfg.get_output_files(
            rules.rank_genes_groups.output.zarr,
            subset_dict=dict(wildcards),
            all_params=True,
        )
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/collect.py'
