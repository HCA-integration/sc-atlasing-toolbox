# HCA Integration Toolbox :toolbox:

**Toolbox of Snakemake pipelines for easy-to-use analyses and benchmarks for building integrated atlases**

- [:rocket: Getting Started](#🚀-getting-started)
- [:gear: Extended Configuration](#⚙️-extended-configuration)
- [:hammer_and_wrench: Trouble Shooting](#🛠️-troubleshooting)

This toolbox provides multiple modules that can be easily combined into custom workflows that leverage the file management of [Snakemake](https://snakemake.readthedocs.io/en/v7.31.1/).
This allows for an efficient and scalable way to run analyses on large datasets that can be easily configured by the user.

The modules are located under `workflow/` and can be run independently or combined into a more complex workflow.
Modules include:

* `load_data`: Loading datasets from URLs and converting them to AnnData objects
* `exploration`, `batch_analysis`, `qc`, `doublets`: Exploration and quality control of datasets
* `merge`, `filter`, `subset`, `relabel`, `split_data`: Basic data manipulation tasks
* `preprocessing`: Preprocessing of datasets (normalization, feature selection, PCA, kNN graph, UMAP)
* `integration`: Running single cell batch correction methods of datasets
* `metrics`: Calculating scIB metrics, mainly for benchmarking of integration methods
* `label_harmonisation`: Provide alignment between unharmonized labels using CellHint


<details>
  <summary>How do I specify a workflow?</summary>

  The heart of the configuration is captured in a YAML (or JSON) configuration file.
  You can find all the modules under `workflow/` and example configuration files under `configs/`.
  Here is an example on a workflow containing the `preprocessing`, `integration` and `metrics` modules:

  ```yaml
  out_dir: /path/to/output/directory
  images: /path/to/image/directory

  DATASETS:

    my_dataset: # custom task/workflow name

      # input specification: map of module name to map of input file name to input file path
      input:
        preprocessing:
          file_1: file_1.h5ad
          file_2: file_2.zarr
        integration: preprocessing # all outputs of module will automatically be used as input
        metrics: integration
      
      # module configuration
      preprocessing:
        highly_variable_genes:
          n_top_genes: 2000
        pca:
          n_pcs: 50
        assemble:
          - normalize
          - highly_variable_genes
          - pca
      
      # module configuration
      integration:
        raw_counts: raw/X
        norm_counts: X
        methods:
          unintegrated:
          scanorama:
            batch_size: 100
          scvi:
            max_epochs: 10
            early_stopping: true

      # module configuration
      metrics:
        unintegrated: layers/norm_counts
        methods:
          - nmi
          - graph_connectivity
  ```

  :sparkling_heart: Beautiful, right? [Read more](#configure-your-workflow) on how configuration works.

</details>


## :rocket: Getting started

### Clone the repository

Depending on whether you have set up SSH or HTTPS with PAT, you can clone the repository with

SSH:
```commandline
git clone git@github.com:HCA-integration/hca_integration_toolbox.git
```

HTTPS:
``` clone
git clone git@github.com:HCA-integration/hca_integration_toolbox.git
```

### Requirements

* Linux (preferred) or MacOS on Intel (not rigorously tested, some bioconda dependencies might not work)
* conda e.g. via [miniforge](https://github.com/conda-forge/miniforge)(recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)


> :memo: **Note** The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [mamba](https://mamba.readthedocs.io).
Please ensure that you have either the mamba or conda pacakage managers installed.
If you use conda, but have never used mamba, consider installing the mamba package into your base environment and use it for all installation commands.
You can still replace all mamba commands with conda commands if you don't install mamba.

### Install conda environments

The different parts of the workflow (modules, rules) require specific conda environments.
The simplest way to install all environments is to run the following script:

```commandline
bash envs/install_all_envs.sh -h # help message for customization
bash envs/install_all_envs.sh
```

> :memo: **Notes**
> 1. The script will create new environments in the `envs` directory if they don't yet exist and update any pre-existing environments.
> 2. If an environment creation fails, the script will skip that environment and you might need to troubleshoot the installation manually.
> 3. The environment names correspond the theire respective file names and are documented under the `name:` directive in the `envs/<env_name>.yaml` file.

If you know you only need certain environments (you can get that information from the README of the module you intend to use), you can install that environment directly.
You will at least require the snakemake environment.

```commandline
mamba env create -f envs/snakemake.yaml
mamba env create -f envs/<env_name>.yaml
```

### Configure your workflow

Configuring your workflow requires setting global variables as well as subworkflows consisting of modules.

#### Configure modules

You can select and combine modules to create a custom workflow by specifying the input and module configuration in a YAML file.
Each instance of a workflow needs a unique task name and it can take any number of inputs consist of modules.

```yaml
DATASETS: # TODO: rename to TASKS

  my_dataset: # custom task/workflow name
    # input specification: map of module name to map of input file name to input file path
    input:
      preprocessing:
        file_1: file_1.h5ad
        file_2: file_2.zarr
      integration: preprocessing # all outputs of module will automatically be used as input
      metrics: integration

  another_dataset:
    ...
 ```

> :warning: **Warning** There can only be one instance of a module as a key in the input mapping (in the backend this is a dictionary). But you can reuse the same module output as input for multiple other modules. The order of the entries in the input mapping doesn't matter. 
s
You can configure the behaviour of each module by specifying their parameters under the same dataset name.
 ```yaml
DATASETS:
  my_dataset:
    input:
      ...

    # module configuration
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
      pca:
        n_pcs: 50
      assemble:
        - normalize
        - highly_variable_genes
        - pca
    
    # module configuration
    integration:
      raw_counts: raw/X
      norm_counts: X
      methods:
        unintegrated:
        scanorama:
          batch_size: 100
        scvi:
          max_epochs: 10
          early_stopping: true

    # module configuration
    metrics:
      unintegrated: layers/norm_counts
      methods:
        - nmi
        - graph_connectivity
```

Each module has a specific set of parameters that can be configured.
Read more about the specific parameters in the README of the module you want to use.

> :memo: **Note** The recommended way to manage your workflow configuration files is to save them outside of the toolbox directory in a directory dedicated to your project. That way you can guarantee the separatation of the toolbox and your own configuration.


#### Global configuration

TODO
* input/output location
* computational resources
* set a Snakemake profile


### Run the pipeline

Before running the pipeline, you need to activate your snakemake environment.

```commandline
conda activate snakemake
```

<details>
  <summary>How does Snakemake work?</summary>

  > The general command for running a pipeline is:
  >
  > ```commandline
  > snakemake <snakemake args>
  > ```
  >
  > The most relevant snakemake arguments are:
  > 
  > + `-n`: dryrun
  > + `--use-conda`: use rule-specific conda environments to ensure all dependencies are met
  > + `-c`: maximum number of cores to be used
  > + `--configfile`: specify a config file to use. The overall workflow already defaults to the config file under `configs/config.yaml`
  > 
  > :bulb: Check out the [snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html) for more commandline arguments.

</details>

#### Create a wrapper script (recommended)
```bash
#!/usr/bin/env bash
set -e -x

snakemake \
  --profile .profiles/czbiohub \
  --configfile \
    configs/computational_resources/czbiohub.yaml \
    configs/integration/config.yaml \
    $@
```

#### Call the pipeline

Dryrun:
```commandline
snakemake -n
```

Run full pipelin with 10 cores:
```commandline
snakemake --use-conda -c10
```

#### Specify which subworkflow you want to run
```commandline
snakemake -l
snakemake load_data_all --use-conda -n
```

## :gear: Extended configuration

### Set defaults
TODO

### Automatic environment management
Snakemake supports automatically creating conda environments for each rule.

```yaml
env_mode: from_yaml
```

### Working with GPUs and Conda environments

Some scripts can run faster if their dependencies are installed with GPU support.
Currently, whether the GPU version of a package with GPU support is installed, depends on the architecture of the system that you install you install the environment on.
That means that if you want to your code to reconnise GPUs when working on a cluster, please make sure you install the conda environments from a node that has access to a GPU.
You can install all missing dependencies in advance:

```
snakemake --use-conda --conda-create-envs-only --cores 1
```

In case you already have an environment installed that doesn't recognise the GPU on a GPU node, gather the environment name from the Snakemake log, remove it manually and then call the pipeline again with `--use-conda` or a corresponding profile.
Snakemake should automatically reinstall that environment.

### Working with CPUs only

If your system doesn't have any GPUs, you can set the following flag in your config:

```yaml
use_gpu: false
```


### Use Snakemake profiles

Different [Snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) are provided
under `.profiles`.
These save defaults for commandline parameters and simplify the snakemake call.
To use a profile e.g. the local profile, call

```commandline
snakemake --profile .profiles/local
```

### Cluster support
TODO

## :hammer_and_wrench: Troubleshooting
TODO