use rule run_method from integration as integration_per_lineage_run with:
    input:
        zarr=rules.preprocessing_per_lineage_assemble.output.zarr,
    output:
        zarr=directory(out_dir / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'adata.zarr'),
        model=touch(directory(out_dir / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'model'))
    benchmark:
        out_dir / 'integration' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        raw_counts=lambda wildcards: get_params(wildcards,parameters,'raw_counts'),
        output_type=lambda wildcards: get_params(wildcards,parameters,'output_type'),
        hyperparams=lambda wildcards: get_params(wildcards,parameters,'hyperparams_dict'),
        env=lambda wildcards: get_params(wildcards,parameters,'env'),
    resources:
        partition=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='partition'),
        qos=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='mem_mb', attempt=attempt),
        gpu=lambda w: get_resource(config,profile=get_params(w,parameters,'resources'),resource_key='gpu'),
        time="2-00:00:00",


use rule postprocess from integration as integration_per_lineage_postprocess with:
    input:
        zarr=rules.integration_per_lineage_run.output.zarr
    output:
        zarr=directory(out_dir / 'postprocess' /paramspace.wildcard_pattern / 'lineage~{lineage}' / 'postprocessed.zarr'),
    resources:
        partition=lambda w: get_resource(config,profile='gpu',resource_key='partition'),
        qos=lambda w: get_resource(config,profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),
