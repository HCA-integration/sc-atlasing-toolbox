# Changelog

## 7.11.2023 Optimise batch PCR analysis

- parallelise permutations per covariate withing script
- use unfiltered files for merging DCP columns
- optimise Preprocessing io

## 19.10.2023 Handle multiple inputs

### New feature: allow multiple inputs to a module
Allow the user to define multiple inputs for a given task ("dataset").

```yaml
DATASETS:
  test_data:
    input:
      preprocessing: 'data.h5ad'
```

The config still work, but it the output files will now have an additional hash code for the input file.

Additionally, the user can now define a mapping of an input file name and a mapping

```yaml
DATASETS:
  test_data:
    input:
      preprocessing:
        file1: 'data.h5ad'
        file2: 'data2.h5ad'
```


This gets particularly handy, when channeling multiple outputs of a module to another module


```yaml
DATASETS:
  test_data:
    input:
      integration: 'data.h5ad'  # provides multiple output files
      metrics: integration  # will take all the output files of integration as input with human readable input file IDs
```


### Breaking changes

1. The output folder structure will look very different now. In case you are using a custom `outputs.yaml` file, please make sure to adapt it to the latest file patterns in `configs/outputs.yaml`

2. The module `integration_per_lineage` is not supported at the moment and will be replaced by a `split_dataset` module combined with `integration`.

### Under the hood

New classes (`workflow/utils/ModuleConfig.py`, `workflow/utils/WildcardParameters.py`, `workflow/utils/InputFiles.py`, `workflow/integration/IntegrationConfig.py`) are now created to handle config inputs and parse wildcard combinations for the workflow.
They make writing new modules easier and improve code reuse.
Most of the functions in `workflow/utils/` will become redundant (they're still included just in case for now).

For example, the integration `Snakefile` goes from this:

<details>
<summary>Old code</summary>

```python

from utils.misc import all_but, unique_dataframe
from utils.config import get_hyperparams, get_resource, get_params_from_config, set_defaults, get_datasets_for_module, get_for_dataset
from utils.wildcards import expand_per, get_params, get_wildcards, wildcards_to_str
from utils.environments import get_env

module_name = 'integration'
config = set_defaults(config,module_name)
out_dir = Path(config['output_dir']) / module_name
image_dir = Path(config['images']) / module_name

# ... 

parameters = pd.read_table(workflow.source_path('params.tsv'))
parameters['output_type'] = parameters['output_type'].str.split(',')
parameters = get_params_from_config(
    config=get_datasets_for_module(config, module_name),
    module_name=module_name,
    config_params=['methods', 'label', 'batch', 'norm_counts', 'raw_counts'],
    wildcard_names=['dataset', 'method', 'label', 'batch', 'norm_counts', 'raw_counts'],
    defaults=config['defaults'],
    explode_by=['method', 'batch'],
).merge(parameters,on='method')

# subset to datasets that have module defined
parameters = parameters[~parameters['method'].isnull()]

# TODO: remove redundant wildcards
# parameters['label'] = np.where(parameters['use_cell_type'], parameters['label'], 'None')
# parameters = unique_dataframe(parameters)

hyperparams_df = get_hyperparams(config,module_name=module_name)
parameters = parameters.merge(hyperparams_df,on=['dataset', 'method'],how='left')
wildcard_names = ['dataset', 'batch', 'label', 'method', 'hyperparams']

# write hyperparameter mapping
Path(out_dir).mkdir(parents=True, exist_ok=True)
unique_dataframe(
    hyperparams_df[['method', 'hyperparams', 'hyperparams_dict']]
).to_csv(out_dir / 'hyperparams.tsv', sep='\t', index=False)

paramspace = Paramspace(
    parameters[wildcard_names],
    filename_params=['method', 'hyperparams'],
    filename_sep='--',
)
```

</details>


to this:

<details>
<summary>New code</summary>

```python
from utils.environments import get_env
from IntegrationConfig import IntegrationConfig

mcfg = IntegrationConfig(
    module_name='integration',
    config=config,
    parameters=workflow.source_path('params.tsv'),
    config_params=['methods', 'batch', 'label', 'norm_counts', 'raw_counts'],
    wildcard_names=['method', 'batch', 'label'],
    rename_config_params={'methods': 'method'},
    explode_by=['method', 'batch'],
)

out_dir = mcfg.out_dir
image_dir = mcfg.image_dir
paramspace = mcfg.get_paramspace()
```
</details>