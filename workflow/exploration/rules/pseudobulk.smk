rule pseudobulk:
    input:
        h5ad=rules.load_data_filter.output.h5ad
    output:
        pca_1_2=out_dir / 'pseudobulk' / '{study}_1_2.png',
        pca_2_3=out_dir / 'pseudobulk' / '{study}_2_3.png',
        pca_scree=out_dir / 'pseudobulk' / '{study}_scree.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        bulk_by='sample',
        color=['donor', 'sex', 'disease'],
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/pseudobulk.py'


use rule pseudobulk as pseudobulk_organ with:
    input:
        h5ad=rules.load_data_merge_organ.output.h5ad
    params:
        dataset=lambda wildcards: wildcards.organ,
        bulk_by='sample',
        color=['dataset', 'reference', 'sex', 'disease', 'assay', 'modalities', 'suspension_type', 'pipeline_version'],
    output:
        pca_1_2=out_dir / 'pseudobulk' / 'organ' / '{organ}_1_2.png',
        pca_2_3=out_dir / 'pseudobulk' / 'organ' / '{organ}_2_3.png',
        pca_scree=out_dir / 'pseudobulk' / 'organ' / '{organ}_scree.png',
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')


rule pseudobulk_all:
    input:
        expand(rules.pseudobulk.output,study=dataset_df['study'].unique()),
        expand(rules.pseudobulk_organ.output,organ=dataset_df['organ'].unique()),
