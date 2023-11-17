rule qc_metrics:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        obs=mcfg.out_dir / params.wildcard_pattern / 'qc_metrics.tsv'
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/qc_metrics.py'


rule qc_metrics_plot:
    input:
        obs=rules.qc_metrics.output.obs
    output:
        joint=mcfg.image_dir / params.wildcard_pattern / 'joint.png',
        joint_mito=mcfg.image_dir / params.wildcard_pattern / 'joint_mito.png',
        joint_log=mcfg.image_dir / params.wildcard_pattern / 'joint_log.png',
        violin=mcfg.image_dir / params.wildcard_pattern / 'violin.png',
        average_jitter=mcfg.image_dir / params.wildcard_pattern / 'average_jitter.png',
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'donor'),
        sample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/qc_metrics_plot.py'


rule doublets:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        obs=mcfg.out_dir / params.wildcard_pattern / 'doublets.tsv',
        zarr=directory(mcfg.out_dir / params.wildcard_pattern / 'doublets.zarr'),
    conda:
        get_env(config, 'qc')
    params:
        batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch'),
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/doublets.py'


rule qc_metrics_all:
    input:
        mcfg.get_output_files(rules.qc_metrics_plot.output),
        # mcfg.get_output_files(rules.doublets.output),
