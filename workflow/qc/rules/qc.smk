rule metrics:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr')
    params:
        gauss_threshold=0.05,
    conda:
        get_env(config, 'qc')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/metrics.py'


rule metrics_all:
    input:
        mcfg.get_output_files(rules.metrics.output),


rule plot_joint:
    input:
        zarr=rules.metrics.output.zarr
    output:
        joint=directory(mcfg.image_dir / params.wildcard_pattern / 'joint_plots'),
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hue', default=[]),
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
    conda:
        get_env(config, 'plots')
    script:
        '../scripts/plot_joint.py'


rule plot_removed:
    input:
        zarr=rules.metrics.output.zarr
    output:
        plots=directory(mcfg.image_dir / params.wildcard_pattern / 'removed'),
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hue', default=[]),
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
    conda:
        get_env(config, 'plots')
    script:
        '../scripts/plot_removed.py'


rule plots_all:
    input:
        mcfg.get_output_files(rules.plot_joint.output),
        mcfg.get_output_files(rules.plot_removed.output),
