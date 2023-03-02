from utils.wildcards import wildcards_to_str


rule barplot:
    input:
        tsv='test/data/integration.benchmark.tsv'
    output:
        png='{out_dir}/barplot_{metric}.png',
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        hue=None,
        facet_row='dataset',
        facet_col=None,
        title='Barplot',
        description=wildcards_to_str,
        dodge=True,
        xlim=(-.01, None),
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/barplot.py'


rule swarmplot:
    input:
        tsv='test/data/metrics.tsv'
    output:
        png='{out_dir}/swarmplot_{metric}.png',
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        title='Swarmplot',
        description=wildcards_to_str,
        ylim=(-.05, 1.05),
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/swarmplot.py'
