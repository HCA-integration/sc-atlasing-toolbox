rule metrics:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        obs=mcfg.out_dir / params.wildcard_pattern / 'qc_metrics.tsv'
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/metrics.py'


rule metrics_plot:
    input:
        obs=rules.metrics.output.obs
    output:
        joint=directory(mcfg.image_dir / params.wildcard_pattern / 'joint_plots'),
        violin=mcfg.image_dir / params.wildcard_pattern / 'violin.png',
        average_jitter=mcfg.image_dir / params.wildcard_pattern / 'average_jitter.png',
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hue', default=mcfg.get_from_parameters(wildcards, 'donor')),
        sample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/metrics_plot.py'


rule metrics_all:
    input:
        mcfg.get_output_files(rules.metrics.output),
        mcfg.get_output_files(rules.metrics_plot.output),
