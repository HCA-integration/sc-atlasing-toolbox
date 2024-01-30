def get_annotated_study(wildcards):
    if wildcards.study in dcp_studies:
        return dict(rules.add_dcp_metadata.output)
    return dict(rules.load_data_merge_study.output)


use rule filter from load_data as load_data_filter_study with:
    input:
        unpack(get_annotated_study)
    output:
        zarr=directory(out_dir / 'filtered' / '{study}.zarr'),
        X=directory(out_dir / 'filtered' / '{study}.zarr' / 'X'),
        removed=directory(out_dir / 'filtered' / 'removed' / '{study}.zarr'),
    params:
        filter=lambda wildcards: config['filter_per_study'][wildcards.study],
        backed=True,
        dask=True,
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,


rule filter_all:
    input: expand(rules.load_data_filter_study.output,**get_wildcards(dataset_df,['study']))
