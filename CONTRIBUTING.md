# Contributing to HCA Integration Toolbox

## Installation

Install the requirements as instructed in the [README](README.md).

## Developing a module

### Module file structure
TODO

### Run a single module

The general command from the directory of interest is:

```commandline
snakemake --configfile <configfile> <snakemake args>
```

You need to specify a config file that is specific to the data you want to run the pipeline on.
This is most useful for testing or reusing the modules for other workflows.

### Testing a module
Make sure a test dataset `pbmc68k.h5ad` exists in `data/`.
This will be the default test object
If it doesn't exist, create it by running `data/generate_data.py` in the `scanpy` environment (see `envs` for environment config files).
Check out the `test` directories of each module for an example.
Here's an example for calling a module pipeline from the repository root:

```commandline
bash workflow/integration/test/run_test.sh
```

### Setting up a new module

## Create a Pull request
...