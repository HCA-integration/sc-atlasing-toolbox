if 'DATASETS' not in config.keys():
    config["DATASETS"] = {}

for study in dataset_df['study']:
    config["DATASETS"][study] = {
        'input': {
            'preprocessing': expand(rules.load_data_filter.output.zarr,study=study)[0],
        },
        'sample': 'sample',
        'label': 'cell_type',
        'batch': 'dataset',
        'lineage' : None,
    }


resource_mode = 'cpu'

rule batch_pcr:
    input:
        zarr=lambda wildcards: expand(rules.preprocessing_pca.output,dataset=wildcards.study)[0],
        metadata=rules.load_data_obs_merge_dcp.output.obs,
    output:
        barplot=out_dir / 'batch_pcr' / '{study}.png',
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
            'library_preparation_protocol.library_construction_method.ontology_label',
            'library_preparation_protocol.end_bias',
            'library_preparation_protocol.strand',
            'library_preparation_protocol.umi_barcode.barcode_read',
            'library_preparation_protocol.umi_barcode.barcode_offset',
            'library_preparation_protocol.umi_barcode.barcode_length',
            'collection_protocol.method.ontology',
            'collection_protocol.method.ontology_label',
            'enrichment_protocol.method.ontology',
            'enrichment_protocol.method.ontology_label',
            'cell_suspension.biomaterial_core.biomaterial_id'
        ],
        permutation_covariates=[
            'assay'
        ],
        sample_key='sample'
    conda:
        '../envs/scib_accel.yaml' if 'os' in config.keys() and config['os'] == 'intel' else '../envs/scib.yaml'
    resources:
        partition=lambda w: get_resource(config,profile=resource_mode,resource_key='partition'),
        mem_mb=get_resource(config,profile=resource_mode,resource_key='mem_mb'),
        qos= lambda w: get_resource(config,profile=resource_mode,resource_key='qos'),
        gpu=lambda w: get_resource(config,profile=resource_mode,resource_key='gpu'),
    script:
        '../scripts/batch_pcr.py'


rule batch_pcr_all:
    input:
        expand(rules.batch_pcr.output,study=dataset_df['study'].unique()),
