rule plot_user:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        dotplot=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'user_markers'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot', default={}),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'layer', default='X'),
        markers=lambda wildcards: get_marker_gene_set(mcfg, wildcards),
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
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
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/plot.py'
