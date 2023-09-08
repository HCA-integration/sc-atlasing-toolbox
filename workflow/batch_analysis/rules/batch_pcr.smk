for dataset in datasets:
    if dataset in config['DATASETS']:
        config["DATASETS"][dataset]['input'].update(
            {
                'preprocessing': get_for_dataset(config, dataset, query=['input', module_name])
            }
        )
        if 'preprocessing' not in config["DATASETS"][dataset]:
            config["DATASETS"][dataset]['preprocessing'] = {}


module preprocessing:
    snakefile: "../../preprocessing/Snakefile"
    config: config

use rule * from preprocessing as preprocessing_*


checkpoint determine_covariates:
    input:
        anndata=lambda w: get_for_dataset(config, w.dataset, ['input', module_name])
    output:
        directory(out_dir / 'batch_pcr' / '{dataset}' / 'covariate_setup')
    params:
        covariates=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'covariates'], default=[]),
        permute_covariates=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'permute_covariates'], default=[]),
        n_permute=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'n_permutations'], default=0),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/determine_covariates.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.determine_covariates.get(**wildcards).output[0]}/{{covariate}}.yaml'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        covariate=glob_wildcards(checkpoint_output).covariate,
        allow_missing=True
    )


rule batch_pcr:
    input:
        anndata=rules.preprocessing_pca.output.zarr,
        setup=get_checkpoint_output,
    output:
        tsv=out_dir / 'batch_pcr' / '{dataset}' / '{covariate}.tsv',
    params:
        n_permute=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'n_permutations'], default=0),
        sample_key=lambda w: get_for_dataset(config, w.dataset, query=[module_name, 'sample_key']),
    conda:
        get_env(config, 'scib_accel')
    resources:
        partition=lambda w: get_resource(config,profile='cpu',resource_key='partition'),
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        qos= lambda w: get_resource(config,profile='cpu',resource_key='qos'),
    script:
        '../scripts/batch_pcr.py'


rule collect:
    input:
        tsv=lambda w: get_from_checkpoint(w, rules.batch_pcr.output.tsv),
    output:
        tsv=out_dir / 'batch_pcr' / '{dataset}.tsv',
    run:
        dfs = [pd.read_table(file) for file in input.tsv]
        df = pd.concat(dfs, ignore_index=True)
        df.to_csv(output.tsv, sep='\t', index=False)


rule plot:
    input:
        tsv=rules.collect.output.tsv,
    output:
        barplot=images_dir / 'batch_pcr' / '{dataset}_bar.png',
        violinplot=images_dir / 'batch_pcr' / '{dataset}_violin.png',
    conda:
        get_env(config, 'plots')
    script:
        '../scripts/plot.py'