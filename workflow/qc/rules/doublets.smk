rule doublets:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        obs=mcfg.out_dir / params.wildcard_pattern / 'doublets.tsv',
        zarr=directory(mcfg.out_dir / 'doublets' / f'{params.wildcard_pattern}.zarr'),
    conda:
        get_env(config, 'qc')
    params:
        batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch'),
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/doublets.py'


rule doublets_all:
    input:
        mcfg.get_output_files(rules.doublets.output),