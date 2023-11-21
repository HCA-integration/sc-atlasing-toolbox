def get_annotated_study(wildcards):
    if wildcards.study in dcp_studies:
        return rules.add_dcp_metadata.output.zarr
    return rules.load_data_merge_study.output.zarr


use rule filter from load_data as load_data_filter_study with:
    input:
        zarr=get_annotated_study
    output:
        zarr=directory(out_dir / 'filtered' / '{study}.zarr'),
        removed=directory(out_dir / 'filtered' / 'removed' / '{study}.zarr'),
    params:
        filter=lambda wildcards: config['filter_per_study'][wildcards.study]
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,


rule filter_all:
    input: expand(rules.load_data_filter_study.output,**get_wildcards(dataset_df,['study']))
