if 'DATASETS' not in config.keys():
    config["DATASETS"] = {}

for study in dataset_df['study']:
    config["DATASETS"][study] = {
        'input': {
            'preprocessing': expand(rules.load_data_filter.output.zarr,study=study)[0],
        },
        'batch': 'dataset',
        'lineage' : None,
    }


resource_mode = 'cpu'

rule batch_pcr:
    input:
        zarr=lambda wildcards: expand(rules.preprocessing_pca.output,dataset=wildcards.study)[0]
    output:
        barplot=out_dir / 'batch_pcr' / '{study}.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        covariates=['sample', 'donor', 'batch', 'assay', 'sex', 'disease', 'self_reported_ethnicity', 'development_stage'],
        permutation_covariates=['assay'],
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
