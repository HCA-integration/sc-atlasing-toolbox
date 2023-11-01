download_columns = ['url', 'dataset', 'collection_id', 'dataset_id', 'project_uuid']

rule download:
    message:
        """
        Downloading file for {wildcards}
        """
    output:
        h5ad=out_dir / 'download' / '{dataset}.h5ad'
    params:
        dataset_df=lambda w: dataset_df.query(f'dataset == "{w.dataset}"')[download_columns].reset_index(drop=True)
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


def get_files_for_metadata_harm(wildcards):
    meta = unlist_dict(
            get_wildcards(dataset_df, columns=all_but(dataset_df.columns,'subset'), wildcards=wildcards)
        )
    file_path = meta['url']
    anno_file = meta['annotation_file']

    # download file if not present
    if not Path(str(file_path)).exists():
        file_path = rules.download.output.h5ad

    # include annotation file if specified
    if pd.notna(anno_file) and pd.notnull(anno_file):
        return {
            'h5ad': file_path,
            'schema': config["schema_file"],
            'annotation_file': anno_file,
        }

    # default files
    return {
        'h5ad': file_path,
        'schema': config["schema_file"],
    }


rule harmonize_metadata:
    input:
        unpack(get_files_for_metadata_harm)
    params:
        meta=lambda wildcards: unlist_dict(
            get_wildcards(dataset_df, columns=all_but(dataset_df.columns,'subset'), wildcards=wildcards)
        )
    output:
        zarr=directory(out_dir / 'harmonize_metadata' / '{dataset}.zarr'),
        plot=image_dir / 'harmonize_metadata' / 'counts_sanity--{dataset}.png',
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,
    # shadow: 'shallow'
    script:
        '../scripts/harmonize_metadata.py'


rule harmonize_metadata_all:
    input:
        expand(rules.harmonize_metadata.output,dataset=dataset_df['dataset'])


use rule merge from load_data as load_data_merge_study with:
    input:
        lambda wildcards: expand(
            rules.harmonize_metadata.output.zarr,
            dataset=get_wildcards(dataset_df,'dataset',wildcards)['dataset']
        )
    output:
        zarr=directory(out_dir / 'merged' / 'study' / '{study}.zarr'),
    params:
        dataset=lambda wildcards: wildcards.study,
        merge_strategy='inner',
        keep_all_columns=True,
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
