config_exploration = config.copy()

if 'DATASETS' not in config_exploration.keys():
    config_exploration["DATASETS"] = {}

for study in dataset_df['study']:
    config_exploration["DATASETS"][study] = dict(
        input=dict(
            preprocessing=expand(rules.load_data_filter.output.zarr,study=study)[0]
        ),
        preprocessing=dict(
            raw_counts='X',
            batch='dataset',
            lineage=None,
        )
    )


module preprocessing:
    snakefile: "../../preprocessing/Snakefile"
    config: config_exploration

use rule * from preprocessing as preprocessing_*

def get_batch_pcr_input(wildcards):
    adata_file = expand(rules.preprocessing_pca.output,dataset=wildcards.study)[0]
    try:
        # try if metadata DCP file exists
        rules.load_data_obs_merge_dcp.input.dcp(wildcards)
        return dict(
            zarr=adata_file,
            metadata=rules.load_data_obs_merge_dcp.output.obs,
        )
    except IndexError as e:
        return dict(zarr=adata_file)


rule batch_pcr:
    input:
        unpack(get_batch_pcr_input)
    output:
        barplot=images_dir / 'batch_pcr' / '{study}.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        covariates=[
            'sample',
            'donor',
            'batch',
            'assay',
            'sex',
            'disease',
            'self_reported_ethnicity',
            'development_stage',
            'sequencing_protocol.method.ontology',
            'library_preparation_protocol.cell_barcode.barcode_read',
            'library_preparation_protocol.cell_barcode.barcode_offset',
            'library_preparation_protocol.cell_barcode.barcode_length',
            'library_preparation_protocol.input_nucleic_acid_molecule.ontology',
            'library_preparation_protocol.library_construction_method.ontology',
            'library_preparation_protocol.end_bias',
            'library_preparation_protocol.strand',
            'library_preparation_protocol.umi_barcode.barcode_read',
            'library_preparation_protocol.umi_barcode.barcode_offset',
            'library_preparation_protocol.umi_barcode.barcode_length',
            'collection_protocol.method.ontology',
            'enrichment_protocol.method.ontology',
            'cell_suspension.biomaterial_core.biomaterial_id'
        ],
        permutation_covariates=[
            'assay',
            # 'sample',
            'batch',
            # 'donor',
            'sequencing_protocol.method.ontology',
            'library_preparation_protocol.cell_barcode.barcode_read',
            'library_preparation_protocol.cell_barcode.barcode_offset',
            'library_preparation_protocol.cell_barcode.barcode_length',
            'library_preparation_protocol.input_nucleic_acid_molecule.ontology',
            'library_preparation_protocol.library_construction_method.ontology',
            'library_preparation_protocol.end_bias',
            'library_preparation_protocol.strand',
            'library_preparation_protocol.umi_barcode.barcode_read',
            'library_preparation_protocol.umi_barcode.barcode_offset',
            'library_preparation_protocol.umi_barcode.barcode_length',
            'collection_protocol.method.ontology',
            'enrichment_protocol.method.ontology',
            'cell_suspension.biomaterial_core.biomaterial_id'
        ],
        n_permute=10,
        sample_key='sample'
    conda:
        get_env(config, 'scib_accel')
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        qos= lambda w: get_resource(config,profile='cpu',resource_key='qos'),
        gpu=lambda w: get_resource(config,profile='cpu',resource_key='gpu'),
    script:
        '../scripts/batch_pcr.py'


rule batch_pcr_all:
    input:
        expand(rules.batch_pcr.output,study=dataset_df['study'].unique()),
