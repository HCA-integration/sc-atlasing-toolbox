rule summary_stats:
    """
    Summarise
    + number of cells
    + number of cells per sample
    + number of cells per donor
    + disease states
    """
    input:
        zarr=rules.load_data_filter.output.zarr
    output:
        tsv=images_dir / 'summary' / 'datasets' / '{study}.tsv',
        sample=images_dir / 'summary' / 'datasets' / '{study}_sample.png',
        donor=images_dir / 'summary' / 'datasets' / '{study}_donor.png',
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/summary_stats.py'


use rule summary_stats as summary_stats_filtered with:
    input:
        zarr=rules.load_data_filter.output.removed
    output:
        tsv=images_dir / 'summary' / 'datasets' / 'filtered' / '{study}.tsv',
        sample=images_dir / 'summary' / 'datasets' / 'filtered' / '{study}_sample.png',
        donor=images_dir / 'summary' / 'datasets' / 'filtered' / '{study}_donor.png',
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')


rule summary_stats_all:
    input:
        tsv=expand(rules.summary_stats.output.tsv,study=dataset_df['study'].unique()),
    output:
        tsv=images_dir / 'summary' / 'all_datasets.tsv',
        aggregate=images_dir / 'summary' / 'all_datasets_aggregated.tsv',
        png=images_dir / 'summary' / 'all_datasets.png',
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot_summary.py'


rule summary_all:
    input:
        rules.summary_stats_all.output,
        expand(rules.summary_stats_filtered.output.tsv,study=dataset_df['study'].unique())
