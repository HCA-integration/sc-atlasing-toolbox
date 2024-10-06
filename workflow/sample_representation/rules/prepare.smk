rule prepare:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / 'dataset~{dataset}' / 'file_id~{file_id}' / 'var_mask~{var_mask}.zarr'),
    params:
        bulk_by=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        aggregate='mean'
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output),
    localrule: True
