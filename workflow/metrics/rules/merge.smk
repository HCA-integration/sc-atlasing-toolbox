rule merge:
    message:
        """
        Merge all metrics for all datasets and methods
        datasets: {params.wildcards[dataset]}
        methods: {params.wildcards[method]}
        """
    input:
        metrics=lambda wildcards: expand_per(rules.run.output,parameters,wildcards,wildcard_names),
        benchmark=lambda wildcards: expand_per(rules.run.benchmark,parameters,wildcards,wildcard_names),
    output:
        tsv=out_dir / 'metrics.tsv',
    params:
        wildcards=get_wildcards(parameters,wildcard_names)
    conda: '../envs/scanpy.yaml'
    resources:
        mem_mb=1000,
        disk_mb=500
    script: '../scripts/merge.py'


use rule merge as merge_per_dataset with:
    message:
        """
        Merge all metrics for {wildcards}
        """
    input:
        metrics=lambda wildcards: expand_per(rules.run.output,parameters,wildcards,all_but(wildcard_names,'dataset')),
        benchmark=lambda wildcards: expand_per(rules.run.benchmark,parameters,wildcards,all_but(wildcard_names,'dataset')),
    output:
        tsv=out_dir / 'per_dataset' / '{dataset}_metrics.tsv',
    params:
        wildcards=lambda wildcards: get_wildcards(parameters,all_but(wildcard_names,'dataset'),wildcards)


use rule merge as merge_per_method with:
    message:
        """
        Merge all metrics for {wildcards}
        """
    input:
        metrics=lambda wildcards: expand_per(rules.run.output,parameters,wildcards,all_but(wildcard_names,'method')),
        benchmark=lambda wildcards: expand_per(rules.run.benchmark,parameters,wildcards,all_but(wildcard_names,'method')),
    output:
        tsv=out_dir / 'per_method' / '{method}.tsv',
    params:
        wildcards=lambda wildcards: get_wildcards(parameters,all_but(wildcard_names,'method'),wildcards)
