rule pseudobulk:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        pca_1_2=mcfg.image_dir / 'pseudobulk' / params.wildcard_pattern / 'pc_1_2.png',
        pca_2_3=mcfg.image_dir / 'pseudobulk' / params.wildcard_pattern / 'pca_2_3.png',
        pca_scree=mcfg.image_dir / 'pseudobulk' / params.wildcard_pattern / 'pca_screeplot.png',
    params:
        dataset=lambda wildcards: wildcards.file_id,
        bulk_by=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca_colors', default=[]),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/pseudobulk.py'


# use rule pseudobulk as pseudobulk_organ with:
#     input:
#         zarr=rules.load_data_merge_organ.output.zarr
#     params:
#         dataset=lambda wildcards: wildcards.organ,
#         bulk_by='sample',
#         color=[
#             'dataset',
#             'study',
#             'reference',
#             'sex',
#             'disease',
#             'assay',
#             'modalities',
#             'pipeline_version',
#             'institution',
#             'self_reported_ethnicity',
#             'development_stage'
#         ],
#     output:
#         pca_1_2=images_dir / 'pseudobulk' / 'organ' / '{organ}_1_2.png',
#         pca_2_3=images_dir / 'pseudobulk' / 'organ' / '{organ}_2_3.png',
#         pca_scree=images_dir / 'pseudobulk' / 'organ' / '{organ}_scree.png',
#     resources:
#         mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')


rule pseudobulk_all:
    input:
        mcfg.get_output_files(rules.pseudobulk.output),