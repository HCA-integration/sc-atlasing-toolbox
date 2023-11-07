rule filter_all:
    input: expand(rules.load_data_filter_study.output,**get_wildcards(dataset_df,['study']))


use rule merge from load_data as load_data_merge_organ with:
    input:
        lambda wildcards: expand(
            rules.load_data_filter_study.output.zarr,
            **get_wildcards(dataset_df,['study'],wildcards),
        ),
    output:
        zarr=directory(out_dir / 'merged' / 'organ' / '{organ}.zarr')
    params:
        dataset=lambda wildcards: wildcards.organ,
        merge_strategy='outer'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    threads:
        dataset_df['dataset'].nunique()


use rule merge from load_data as load_data_merge_organ_filter with:
    input:
        lambda wildcards: expand(
            rules.load_data_filter_study.output.removed,
            **get_wildcards(dataset_df,['study'],wildcards),
        ),
    output:
        zarr=directory(out_dir / 'merged' / 'organ' / 'filtered' / '{organ}.zarr')
    params:
        dataset=lambda wildcards: wildcards.organ,
        merge_strategy='outer'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    threads:
        dataset_df['dataset'].nunique()


use rule merge from load_data as load_data_merge_subset with:
    input:
        lambda wildcards: expand(
            rules.load_data_filter_study.output.zarr,
            **get_wildcards(dataset_df,['study'],wildcards),
        ),
    output:
        zarr=directory(out_dir / 'merged' / 'subset' / '{organ}-{subset}.zarr')
    params:
        dataset=lambda wildcards: f'{wildcards.organ}-{wildcards.subset}',
        merge_strategy='outer'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    threads:
        dataset_df['dataset'].nunique()


rule  merge_subset_all:
    input:
        expand(
            rules.load_data_merge_subset.output,
            **get_wildcards(dataset_df,['organ', 'subset'],drop_na=True)
        )
