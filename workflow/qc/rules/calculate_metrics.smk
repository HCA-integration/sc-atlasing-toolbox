rule autoqc:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / 'autoqc' / f'{params.wildcard_pattern}.zarr')
    params:
        gauss_threshold=0.05,
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    conda:
        get_env(config, 'qc')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/autoqc.py'


rule autoqc_all:
    input:
        mcfg.get_output_files(rules.autoqc.output),
    localrule: True


rule get_thresholds:
    input: 
        zarr=rules.autoqc.output.zarr
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
        tsv=mcfg.image_dir / params.wildcard_pattern / 'thresholds.tsv',
        qc_stats=mcfg.image_dir / params.wildcard_pattern / 'qc_stats.tsv'
    params:
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
        alternative_thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'alternative_thresholds', default={}),
    conda:
        get_env(config, 'qc')
    localrule: True
    script:
        '../scripts/get_thresholds.py'


rule merge_thresholds:
    input:
        tsv=lambda wildcards: mcfg.get_output_files(
            rules.get_thresholds.output.tsv,
            subset_dict=wildcards
        ),
        qc_stats=lambda wildcards: mcfg.get_output_files(
            rules.get_thresholds.output.qc_stats,
            subset_dict=wildcards
        ),
    output:
        tsv=mcfg.image_dir / 'dataset~{dataset}' / 'thresholds.tsv',
        qc_stats=mcfg.image_dir / 'dataset~{dataset}' / 'qc_stats.tsv',
    localrule: True
    run:
        import pandas as pd
        
        df = pd.concat([pd.read_table(file) for file in input.tsv])
        print(df)
        df.to_csv(output.tsv, sep='\t', index=False)

        df = pd.concat([pd.read_table(file) for file in input.qc_stats])
        print(df)
        df.to_csv(output.qc_stats, sep='\t', index=False)


rule thresholds_all:
    input:
        mcfg.get_output_files(rules.get_thresholds.output),
        mcfg.get_output_files(rules.merge_thresholds.output),
    localrule: True
