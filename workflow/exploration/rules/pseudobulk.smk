rule pseudobulk:
    input:
        zarr=rules.load_data_filter.output.zarr
    output:
        pca_1_2=images_dir / 'pseudobulk' / '{study}_1_2.png',
        pca_2_3=images_dir / 'pseudobulk' / '{study}_2_3.png',
        pca_scree=images_dir / 'pseudobulk' / '{study}_scree.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        bulk_by='sample',
        color=['donor', 'assay', 'sex', 'disease', 'self_reported_ethnicity', 'development_stage', 'batch'],
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/pseudobulk.py'


use rule pseudobulk as pseudobulk_organ with:
    input:
        zarr=rules.load_data_merge_organ.output.zarr
    params:
        dataset=lambda wildcards: wildcards.organ,
        bulk_by='sample',
        color=[
            'dataset',
            'study',
            'reference',
            'sex',
            'disease',
            'assay',
            'modalities',
            'pipeline_version',
            'institution',
            'self_reported_ethnicity',
            'development_stage'
        ],
    output:
        pca_1_2=images_dir / 'pseudobulk' / 'organ' / '{organ}_1_2.png',
        pca_2_3=images_dir / 'pseudobulk' / 'organ' / '{organ}_2_3.png',
        pca_scree=images_dir / 'pseudobulk' / 'organ' / '{organ}_scree.png',
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')


rule pseudobulk_all:
    input:
        expand(rules.pseudobulk.output,study=dataset_df['study'].unique()),
        expand(rules.pseudobulk_organ.output,organ=dataset_df['organ'].unique()),
