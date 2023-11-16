rule qc_stats:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        obs=mcfg.out_dir / params.wildcard_pattern / 'obs.tsv'
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/qc_stats.py'


rule qc_plot:
    input:
        obs=rules.qc_stats.output.obs
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
        '../scripts/qc_plot.py'


rule qc_all:
    input:
        mcfg.get_output_files(rules.qc_plot.output),
