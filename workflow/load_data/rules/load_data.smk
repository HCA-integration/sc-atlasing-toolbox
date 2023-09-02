rule download:
    message:
        """
        Downloading file {wildcards}
        """
    output:
        h5ad=out_dir / 'download' / '{dataset}.h5ad'
    params:
        dataset_df=lambda w: dataset_df.query(f'dataset == "{w.dataset}"').reset_index(drop=True)
    conda:
        get_env(config, 'scanpy', env_dir='../../envs/')
    script:
        'scripts/download.py'

datasets_to_download = dataset_df[
    dataset_df['url'].apply(lambda x: not Path(x).is_file())
]['dataset']
rule download_all:
    input:
        expand(rules.download.output, dataset=datasets_to_download)


def get_h5ad(wildcards):
    file_path = dataset_df[dataset_df['dataset'] == wildcards.dataset].astype(str).reset_index()['url'][0]
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
        get_env(config, 'scanpy', env_dir='../../envs/')
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),
        disk_mb=20000,
    shadow: 'shallow'
    script:
        'scripts/metadata.py'


rule metadata_all:
    input:
        expand(rules.metadata.output,dataset=dataset_df['dataset'])
