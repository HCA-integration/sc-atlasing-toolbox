rule get_thresholds:
    input: 
        zarr=rules.metrics.output.zarr
    output:
        tsv=mcfg.image_dir / params.wildcard_pattern / 'thresholds.tsv'
    params:
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
    conda:
        get_env(config, 'qc')
    script:
        '../scripts/get_thresholds.py'


rule merge_thresholds:
    input:
        tsv=lambda wildcards: mcfg.get_output_files(rules.get_thresholds.output, subset_dict=wildcards),
    output:
        tsv=mcfg.image_dir / 'dataset~{dataset}' / 'thresholds.tsv'
    run:
        import pandas as pd
        
        df = pd.concat([pd.read_table(file) for file in input.tsv])
        print(df)
        df.to_csv(output.tsv, sep='\t', index=False)


rule thresholds_all:
    input: 
        mcfg.get_output_files(rules.get_thresholds.output),
        mcfg.get_output_files(rules.merge_thresholds.output),