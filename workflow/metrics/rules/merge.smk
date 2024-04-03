metric_wildcards = mcfg.get_wildcard_names() + ['metric', 'overwrite_file_id']

rule merge_metrics:
    message:
        """
        Merge all metrics for all datasets and methods
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=wildcards, all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=wildcards, all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'extra_columns.txt',
    params:
        wildcards=mcfg.get_wildcards(as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=mcfg.get_wildcards(as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    conda:
        get_env(config, 'scanpy')
    localrule: True
    resources:
        mem_mb=1000,
        disk_mb=500
    script: '../scripts/merge.py'


use rule merge_metrics as merge_metrics_per_dataset with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_dataset' / '{dataset}_metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_dataset' / '{dataset}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),


use rule merge_metrics as merge_metrics_per_batch with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_batch' / '{batch}_metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_batch' / '{batch}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='batch', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='batch', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),


use rule merge_metrics as merge_metrics_per_label with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_label' / '{label}_metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_label' / '{label}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='label', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='label', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),


use rule merge_metrics as merge_metrics_per_file with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'per_file' / '{file_id}.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'per_file' / '{file_id}_extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='file_id', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='file_id', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),


rule merge_metrics_all:
    input:
        rules.merge_metrics.output,
        mcfg.get_output_files(rules.merge_metrics_per_dataset.output, wildcard_names=['dataset']),
        mcfg.get_output_files(rules.merge_metrics_per_batch.output, wildcard_names=['batch']),
        mcfg.get_output_files(rules.merge_metrics_per_label.output, wildcard_names=['label']),
        mcfg.get_output_files(rules.merge_metrics_per_file.output, wildcard_names=['file_id']),
    localrule: True