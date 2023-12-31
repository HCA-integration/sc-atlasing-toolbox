use rule run_method from integration as integration_run_method with:
    message:
       """
       Integration: Run {wildcards.method} on {wildcards.dataset}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       """
    input:
        h5ad=lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern.replace('--output_type~{output_type}', '') / 'adata.zarr'),
        model=touch(directory(out_dir / paramspace.wildcard_pattern.replace('--output_type~{output_type}', '') / 'model'))
    benchmark:
        out_dir / paramspace.wildcard_pattern.replace('--output_type~{output_type}', '') / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', exclude=['output_type']),
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', exclude=['output_type']),
        output_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_types', exclude=['output_type']),
        hyperparams=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hyperparams_dict', exclude=['output_type']),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env', exclude=['output_type']),
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w)),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w)),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w)),
        time="2-00:00:00",


def update_neighbors_args(wildcards):
    args = mcfg.get_for_dataset(wildcards.dataset, ['preprocessing', 'neighbors'], default={})
    output_type = wildcards.output_type
    if output_type == 'full':
        args |= {'use_rep': 'X_pca'}
    elif output_type == 'embed':
        args |= {'use_rep': 'X_emb'}
    elif output_type == 'knn':
        args = False
    else:
        raise ValueError(f'Unknown output_type {output_type}')
    return args


use rule neighbors from preprocessing as integration_postprocess with:
    input:
        zarr=rules.integration_run_method.output.zarr,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'postprocessed.zarr'),
    params:
        args=update_neighbors_args,
        extra_uns=lambda wildcards: {'output_type': wildcards.output_type},
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),


rule run_all:
    input:
        mcfg.get_output_files(rules.integration_postprocess.output)
