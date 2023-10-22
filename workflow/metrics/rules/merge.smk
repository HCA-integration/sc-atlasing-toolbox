rule merge:
    message:
        """
        Merge all metrics for all datasets and methods
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=wildcards),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=wildcards),
    output:
        tsv=mcfg.out_dir / 'results' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'extra_columns.txt',
    params:
        wildcards=mcfg.get_wildcards(as_df=True),
        wildcards_string=mcfg.get_wildcards(as_df=True).to_string(index=False),
    conda:
        get_env(config, 'scanpy')
    group:
        'metrics_merge'
    resources:
        mem_mb=1000,
        disk_mb=500
    script: '../scripts/merge.py'


use rule merge as merge_per_dataset with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards)),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards)),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_dataset' / '{dataset}_metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_dataset' / '{dataset}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True).to_string(index=False)


use rule merge as merge_per_file with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards)),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards)),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_file' / '{file_id}.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_file' / '{file_id}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='file_id', as_df=True),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='file_id', as_df=True).to_string(index=False)


rule merge_all:
    input:
        mcfg.get_output_files(rules.merge_per_dataset.output, wildcard_names=['dataset']),
        mcfg.get_output_files(rules.merge_per_file.output, wildcard_names=['file_id']),