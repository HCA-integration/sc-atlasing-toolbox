rule qc:
    input:
        zarr=rules.load_data_filter.output.zarr
    output:
        joint=images_dir / 'qc' / '{study}' / 'joint.png',
        joint_log=images_dir / 'qc' / '{study}' / 'joint_log.png',
        violin=images_dir / 'qc' / '{study}' / 'violin.png',
        average_jitter=images_dir / 'qc' / '{study}' / 'average_jitter.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        hue='donor'
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/qc_plot.py'


use rule qc as qc_organ with:
    input:
        zarr=rules.load_data_merge_organ.output.zarr
    output:
        joint=images_dir / 'qc' / 'organ' / '{organ}' / 'joint.png',
        joint_log=images_dir / 'qc' / 'organ' / '{organ}' / 'joint_log.png',
        violin=images_dir / 'qc' / 'organ' / '{organ}' / 'violin.png',
        average_jitter=images_dir / 'qc' / 'organ' / '{organ}' / 'average_jitter.png',
    params:
        dataset=lambda wildcards: wildcards.organ,
        hue='study'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')


use rule qc as qc_filtered with:
    input:
        zarr=rules.merge_organ_filter.output.zarr
    output:
        joint=images_dir / 'qc' / 'organ' / '{organ}' / 'filtered_joint.png',
        joint_log=images_dir / 'qc' / 'organ' / '{organ}' / 'filtered_joint_log.png',
        violin=images_dir / 'qc' / 'organ' / '{organ}' / 'filtered_violin.png',
        average_jitter=images_dir / 'qc' / 'organ' / '{organ}' / 'filtered_average_jitter.png',
    params:
        dataset=lambda wildcards: wildcards.organ,
        hue='study'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')


rule qc_all:
    input:
        expand(rules.qc.output,study=dataset_df['study'].unique()),
        expand(rules.qc_organ.output,organ=dataset_df['organ'].unique()),
        expand(rules.qc_filtered.output,organ=dataset_df['organ'].unique())
