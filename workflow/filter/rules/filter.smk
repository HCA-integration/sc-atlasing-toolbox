rule filter:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
    params:
        remove_by_column=lambda wildcards: mcfg.get_from_parameters(wildcards, 'remove_by_column', default={}),
        backed=lambda wildcards: mcfg.get_from_parameters(wildcards, 'backed'),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask'),
        subset=lambda wildcards: mcfg.get_from_parameters(wildcards, 'subset'),
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt),
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    script:
        '../filter.py'


rule filter_all:
    input: mcfg.get_output_files(rules.filter.output)
