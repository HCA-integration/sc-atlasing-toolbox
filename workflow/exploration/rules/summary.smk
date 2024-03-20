rule summary_stats:
    """
    Summarise
    + number of cells
    + number of cells per sample
    + number of cells per donor
    + disease states
    """
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        tsv=mcfg.image_dir / 'summary' / params.wildcard_pattern / 'stats.tsv',
        sample=mcfg.image_dir / 'summary' / params.wildcard_pattern / 'per_sample.png',
        donor=mcfg.image_dir / 'summary' / params.wildcard_pattern / 'per_donor.png',
    params:
        sample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
        donor=lambda wildcards: mcfg.get_from_parameters(wildcards, 'donor'),
        categories=lambda wildcards: mcfg.get_from_parameters(wildcards, 'summary_columns', default=[]),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/summary_stats.py'


per_dataset_pattern = params.wildcard_pattern.replace('/file_id~{file_id}', '')

rule summary_stats_all:
    input:
        tsv=mcfg.get_output_files(rules.summary_stats.output.tsv, allow_missing=True),
    output:
        tsv=mcfg.image_dir / 'summary' / per_dataset_pattern / 'per_dataset_stats.tsv',
        aggregate=mcfg.image_dir / 'summary' / per_dataset_pattern / 'per_dataset_stats_aggregated.tsv',
        png=mcfg.image_dir / 'summary' / per_dataset_pattern / 'per_dataset_stats.png',
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/summary_plot.py'


rule summary_all:
    input:
        mcfg.get_output_files(rules.summary_stats_all.output, exclude=['file_id']),
