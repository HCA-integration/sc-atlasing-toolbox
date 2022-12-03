# Commonly used workflows

Module for workflows that are commonly used in other modules.

* dependency graph (rule and job graphs)
* adding annotations to obs
* plotting results from TSV file

## Dependency graph
Rules in `rules/dependency_graph.smk`

* `rulegraph`: plot and rule graph for a given target
* `dag` job graph for a given target

No input required, as the graphs are built on the given target.
The graphs are visualised via `dot` and saved as PNG.

## Plots

Rules in `rules/plots.smk`

Facetted plots for comparing computational and quality metrics.
Currently supported plots are:

* `barplot`

### Input

TSV file containing the columns needed by facetting and metrics columns
Which columns are required for facetting the plot can be specified in the `params` directive of the rule.
The rule can be imported as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules)
so that the inputs, outputs and parameters can be modified.

#### Example

The rules expect the benchmark output provided by the Snakemake directive.
The snakemake benchmark columns are as follows:

```tsv
method    s    h:m:s    max_rss    max_vms max_uss    max_pss io_in    io_out    mean_load    cpu_time
```

The entries are explained nicely
here: https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake

An example configuration of the plot parameters:

```snakemake
params:
    metric='s',
    category='method',
    hue=None,
    facet_row=None,
    facet_col=None,
    title='Barplot',
    description=lambda wildcards: ' '.join([f'{key}={value}' for key, value in wildcards.items()]),
    dodge=True,
```

### Output
PNG image of the plot.