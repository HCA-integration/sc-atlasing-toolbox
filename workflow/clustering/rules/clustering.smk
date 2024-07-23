def get_neighbors_file(wildcards):
    level = int(wildcards.get('level', 1))
    if level > 1:
        return expand(
            rules.clustering_cluster.output.zarr,
            level=level - 1,
            allow_missing=True
        )[0]
    if mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return expand(
            rules.clustering_compute_neighbors.output.zarr,
            algorithm='None',
            resolution='None',
            level=1,
            allow_missing=True,
        )[0]
    return mcfg.get_input_file(**wildcards)


use rule neighbors from preprocessing as clustering_compute_neighbors with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'neighbors' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(
            {k: v for k, v in wildcards.items() if k not in ['algorithm', 'resolution']},
            'neighbors',
            default={},
            exclude=['algorithm', 'resolution'],
        ),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'resolutions' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
        neighbors_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
        clustering_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'kwargs', default={}),
    threads:
        lambda wildcards: 5 * int(wildcards.level) - 4
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
