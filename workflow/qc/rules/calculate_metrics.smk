rule metrics:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr')
    params:
        gauss_threshold=0.05,
    conda:
        get_env(config, 'qc')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/metrics.py'


rule metrics_all:
    input:
        mcfg.get_output_files(rules.metrics.output),
