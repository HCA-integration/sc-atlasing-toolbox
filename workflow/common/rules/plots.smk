rule barplot:
    input:
        tsv='test/data/integration.benchmark.tsv'
    output:
        png='{out_dir}/{metric}.png',
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        hue=None,
        facet_row='dataset',
        facet_col=None,
        title='Barplot',
        description=lambda wildcards: ' '.join([f'{key}={value}' for key, value in wildcards.items()]),
        dodge=True,
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/barplot.py'


rule test_plots:
    input:
        expand(rules.barplot.output,out_dir='test/out/benchmark',metric='s'),
        expand(rules.barplot.output,out_dir='test/out/benchmark',metric='max_uss'),
        expand(rules.barplot.output,out_dir='test/out/benchmark',metric='mean_load'),
