rule download:
    message:
        """
        Downloading file for {wildcards}
        """
    output:
        h5ad=out_dir / 'download' / '{dataset}.h5ad'
    params:
        dataset_df=lambda w: dataset_df.query(f'dataset == "{w.dataset}"').reset_index(drop=True)
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    script:
        '../scripts/download.py'


rule download_all:
    input:
        expand(
            rules.download.output,
            dataset=dataset_df[dataset_df['url'].apply(lambda x: not Path(x).is_file())]['dataset']
        )


def get_h5ad(wildcards):
    df = dataset_df.query(f'dataset == "{wildcards.dataset}"').reset_index()
    file_path = df.astype(str)['url'][0]
    if Path(str(file_path)).exists():
        return file_path
    return rules.download.output.h5ad


rule metadata:
    input:
        h5ad=get_h5ad,
        schema=config["schema_file"]
    params:
        meta=lambda wildcards: unlist_dict(
            get_wildcards(dataset_df,columns=all_but(dataset_df.columns,'subset'),wildcards=wildcards)
        )
    output:
        zarr=directory(out_dir / 'processed/' / '{dataset}.zarr'),
        plot=out_dir / 'processed/counts/{dataset}.png',
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,
    # shadow: 'shallow'
    script:
        '../scripts/metadata.py'


rule metadata_all:
    input:
        expand(rules.metadata.output,dataset=dataset_df['dataset'])


use rule merge from load_data as load_data_merge_study with:
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


rule merge_study_all:
    input: expand(rules.load_data_merge_study.output,**get_wildcards(dataset_df,['study']))


use rule filter from load_data as load_data_filter_study with:
    input:
        zarr=rules.load_data_merge_study.output.zarr
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
