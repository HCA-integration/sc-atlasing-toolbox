rule summary_stats:
    """
    Summarise
    + number of cells
    + number of cells per sample
    + number of cells per donor
    + disease states
    """
    input:
        h5ad=rules.load_data_filter.output.h5ad
    output:
        tsv=out_dir / 'summary' / 'datasets' / '{study}.tsv',
        sample=out_dir / 'summary' / 'datasets' / '{study}_sample.png',
        donor=out_dir / 'summary' / 'datasets' / '{study}_donor.png',
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/summary_stats.py'


rule summary_stats_all:
    input:
        tsv=expand(rules.summary_stats.output.tsv,study=dataset_df['study'].unique()),
    output:
        tsv=out_dir / 'summary' / 'all_datasets.tsv',
        aggregate=out_dir / 'summary' / 'all_datasets_aggregated.tsv',
        png=out_dir / 'summary' / 'all_datasets.png',
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/plot_summary.py'
