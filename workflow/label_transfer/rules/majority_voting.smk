rule majority_voting:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
        plots=directory(mcfg.image_dir / f'{paramspace.wildcard_pattern}'),
    params:
        majority_reference=lambda wildcards: mcfg.get_from_parameters(wildcards, 'majority_reference', default={}),
        majority_consensus=lambda wildcards: mcfg.get_from_parameters(wildcards, 'majority_consensus', default={}),
    conda:
        get_env(config, 'scanpy', env_dir='envs')
    resources:
        partition=lambda w: mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w: mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/majority_voting.py'