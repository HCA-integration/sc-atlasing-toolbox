# HCA pipelines

Pipelines for standard analyses of HCA datasets.

## Modules

There are multiple standalone modules that can be run independently or combined.
The modules are located under `workflow`.
The `configs` directory contains files that specify inputs, output location and dataset specific parameters.

+ `config.yaml`: datasets specific parameters to be used for different modules
+ `datsets.tsv`: annotation of datasets of interest and URL. Used for downloading datasets
+ `modules.tsv`: specification of which modules and submodules are to be run

## Setting up the pipeline

Create a conda environment for snakemake

```commandline
mamba env create -f envs/snakemake.yaml
```

> Note: [mamba](https://mamba.readthedocs.io/en/latest/installation.html_) improves install times drastically.
> All mamba commands can be replaced by `conda`.

## Running the pipeline

There are generally two ways to run the pipeline, either from within a module directory (for module-specific runs) or at
the top level (for full pipeline run).

The general command for running a pipeline is:

```commandline
snakemake <snakemake args>
```

The most relevant snakemake arguments are:

+ `-n`: dryrun
+ `--use-conda`: use rule-specific conda environments to ensure all dependencies are met
+ `-c`: maximum number of cores to be used
+ `--configfile`: specify a config file to use. The overall workflow already defaults to the config file
  under `configs/config.yaml`

Refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more
commandline arguments.

### Run the full pipeline

First make sure that all configuration files (under `configs/`) define which datasets you want to run and which modules
are included.
Which modules are computed exactly can be defined in `configs/modules.tsv`.
Then use the following commands to call the pipeline.

Dryrun:
```commandline
snakemake -n
```

Run full pipelin with 10 cores:
```commandline
snakemake --use-conda -c10
```

> Note: the full pipeline doesn't yet combine input and outputs of all modules. That is still WIP


### Run parts of the full pipeline

Check out the `workflow/Snakefile` for the exact parts of the pipeline you want to run.
For example, if you just want the data loading workflow, the rule is called `load_data_all` and the snakemake command
will be:

```commandline
snakemake load_data_all --use-conda -n
```

### Run a single module

The general command from the directory of interest is:

```commandline
snakemake --configfile <configfile> <snakemake args>
```

You need to specify a config file that is specific to the data you want to run the pipeline on.
This is most useful for testing or reusing the modules for other workflows.
Check out the `test` directories of each module for an example.

### Use Snakemake profiles

Different [Snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) are provided
under `.profiles`.
These save defaults for commandline parameters and simplify the snakemake call.
To use a profile e.g. the local profile, call

```commandline
snakemake --profile .profiles/local
```