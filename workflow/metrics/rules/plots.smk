module plots:
    snakefile: "../../common/rules/plots.smk"
    config: config


# barplots

use rule barplot from plots as metrics_barplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / '{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_row='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,


use rule barplot from plots as metrics_barplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=out_dir / 'plots' / 'per_dataset' / '{metric}-barplot-{dataset}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,


use rule barplot from plots as metrics_barplot_per_method with:
    input:
        tsv=rules.merge_per_method.output.tsv
    output:
        png=out_dir / 'plots' / 'per_method' / '{metric}-barplot-{method}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,


# swarm plots

use rule swarmplot from plots as metrics_swarmplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_col='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,


use rule swarmplot from plots as metrics_swarmplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=out_dir / 'plots' / 'per_dataset' / '{metric}-swarmplot-{dataset}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,


use rule swarmplot from plots as metrics_swarmplot_per_method with:
    input:
        tsv=rules.merge_per_method.output.tsv
    output:
        png=out_dir / 'plots' / 'per_method' / '{metric}-swarmplot-{method}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,


rule compare_metrics:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / 'comparison.png'
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/plots/comparison.py'


# Funkyheatmap

rule funkyheatmap:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / 'funky_heatmap.png'
    conda:
        '../envs/plots.yaml'
    singularity:
        'ghcr.io/dynverse/funky_heatmap:latest'
    script:
        '../scripts/plots/funkyheatmap.R'


rule funkyheatmap_standalone:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=out_dir / 'plots' / 'funky_heatmap.png'
    shell:
        """
        wget https://github.com/dynverse/funkyheatmap/releases/latest/download/executable.zip
        unzip -o executable.zip -d funky_heatmap
        ./funky_heatmap/funky_heatmap --data {input.tsv} --output {output.png}
        """
