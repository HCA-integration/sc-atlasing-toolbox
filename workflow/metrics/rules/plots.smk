module plots:
   snakefile: "../../common/rules/plots.smk"
   config: config


def get_merged_tsv(wildcards):
    merge = wildcards.out_dir
    if merge == 'all':
        return rules.merge.output.tsv
    elif merge == 'per_dataset':
        return expand(rules.merge_per_dataset.output,zip,**get_wildcards(parameters,['dataset']))
    elif merge == 'per_method':
        return expand(rules.merge_per_method.output,zip,**get_wildcards(parameters,['method']))
    else:
        raise ValueError(f'invalid merge strategy "{merge}"')


# TODO: deal with per dataset and per method outputs
use rule barplot from plots as metrics_computation_plot with:
    input:
        tsv=get_merged_tsv
    output:
        png=out_dir / 'plots' / '{out_dir}' / 'computation_{metric}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='output_type',
        facet_col='dataset',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,


use rule swarmplot from plots as metrics_swarmplot with:
    input:
        tsv=get_merged_tsv
    output:
        png=out_dir / 'plots' / '{out_dir}' / 'swarmplot_{metric}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='output_type',
        facet_row='dataset',
        title='Metrics',
        description=wildcards_to_str,


rule compare_metrics:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / 'comparison.png'
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/plots/comparison.py'
