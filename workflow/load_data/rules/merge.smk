rule merge_study:
    """
    Merge datasets that had to be split
    """
    input:
        lambda wildcards: expand(
            rules.metadata.output.zarr,
            dataset=get_wildcards(dataset_df,'dataset',wildcards)['dataset']
        )
    output:
        zarr=directory(out_dir / 'merged' / 'study' / '{study}.zarr'),
    params:
        dataset=lambda wildcards: wildcards.study,
        merge_strategy='inner'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    # shadow: 'minimal'
    script:
        'scripts/merge.py'


rule merge_study_all:
    input: expand(rules.merge_study.output,**get_wildcards(dataset_df,['study']))


rule filter:
    """
    QC filter per study
    TODO: keep filtered out cells for stats
    """
    input:
        zarr=rules.merge_study.output.zarr
    output:
        zarr=directory(out_dir / 'filtered' / '{study}.zarr'),
        removed=directory(out_dir / 'filtered' / 'removed' / '{study}.zarr'),
    params:
        filter=lambda wildcards: config['filter_per_study'][wildcards.study]
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    # shadow: 'minimal'
    script:
        'scripts/filter.py'


rule filter_all:
    input: expand(rules.filter.output,**get_wildcards(dataset_df,['study']))


rule merge_organ:
    input:
        lambda wildcards: expand(
            rules.filter.output.zarr,
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
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    # shadow: 'minimal'
    script:
        'scripts/merge.py'


rule merge_organ_filter:
    input:
        lambda wildcards: expand(
            rules.filter.output.removed,
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
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    # shadow: 'minimal'
    script:
        'scripts/merge.py'


rule merge_subset:
    input:
        lambda wildcards: expand(
            rules.filter.output.zarr,
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
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    # shadow: 'minimal'
    script:
        'scripts/merge.py'


rule  merge_subset_all:
    input:
        expand(
            rules.merge_subset.output,
            **get_wildcards(dataset_df,['organ', 'subset'],drop_na=True)
        )
