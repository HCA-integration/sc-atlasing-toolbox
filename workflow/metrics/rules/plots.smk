module common:
    snakefile: "../../common/Snakefile"
    config: config


use rule * from common as common_ *


def get_merged_tsv(wildcards):
    merge = wildcards.merge
    if merge == '':
        return rules.merge.output.tsv
    elif merge == 'per_dataset':
        return rules.merge_per_dataset.output.tsv
    elif merge == 'per_method':
        return rules.merge_per_method.output.tsv
    else:
        raise ValueError(f'invalid merge strategy "{merge}"')


use rule barplot from common as computation_plot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / '{out_dir}computation_{metric}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='output_type',
        facet_col='dataset',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,


use rule swarmplot from common as swarmplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / '{out_dir}swarmplot_{metric}.png'
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
