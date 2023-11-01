# for _, row in mcfg.get_wildcards(as_df=True).iterrows():
#     dataset = row['dataset']
#     file_id = row['file_id']
#     config["DATASETS"][dataset]['input'].update(
#             {
#                 'preprocessing': mcfg.get_input_file(dataset, file_id)
#             }
#     )
#     pp_config = mcfg.get_defaults(module_name='preprocessing')
#     pp_config.update(config["DATASETS"][dataset].get("preprocessing", {}))
#     pp_config.update(dict(raw_counts='X', batch='study'))
#     config["DATASETS"][dataset]["preprocessing"] = pp_config


# module preprocessing:
#     snakefile: "../../preprocessing/Snakefile"
#     config: config

# use rule * from preprocessing as preprocessing_*


checkpoint determine_covariates:
    input:
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        directory(mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr' / 'covariate_setup')
    params:
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'covariates', default=[]),
        permute_covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'permute_covariates', default=[]),
        n_permute=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=10),
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
        anndata=lambda wildcards: mcfg.get_input_file(**wildcards),
        # anndata=rules.preprocessing_pca.output.zarr,
        setup=get_checkpoint_output,
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr' / '{covariate}.tsv',
    params:
        n_permute=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=10, check_query_keys=False),
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key', check_query_keys=False),
    conda:
        get_env(config, 'scib_accel')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=10, check_query_keys=False, as_type=float)
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
    script:
        '../scripts/batch_pcr.py'


rule collect:
    input:
        tsv=lambda wildcards: get_from_checkpoint(wildcards, rules.batch_pcr.output.tsv),
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr.tsv',
    run:
        dfs = [pd.read_table(file) for file in input.tsv]
        if len(dfs) == 0:
            df = pd.DataFrame(columns=['covariate', 'pcr', 'permuted', 'n_covariates'])
        else:
            df = pd.concat(dfs, ignore_index=True)
        df.to_csv(output.tsv, sep='\t', index=False)


rule plot:
    input:
        tsv=rules.collect.output.tsv,
    output:
        barplot=mcfg.image_dir / paramspace.wildcard_pattern / 'batch_pcr_bar.png',
        violinplot=mcfg.image_dir / paramspace.wildcard_pattern / 'batch_pcr_violin.png',
    conda:
        get_env(config, 'plots')
    script:
        '../scripts/plot.py'