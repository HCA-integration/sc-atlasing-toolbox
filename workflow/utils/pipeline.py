"""
Utils that are specific to the overall pipeline but not the modules
"""

from pathlib import Path

from .config import _get_or_default_from_config, set_defaults
from .ModuleConfig import ModuleConfig
from .InputFiles import InputFiles


def update_module_configs(config, params):
    """
    Update config with parameters from modules TSV
    """
    # Process params per module
    for dataset in config['DATASETS'].keys():

        def _get(params, dataset, module):
            submodules = params.query('dataset == @dataset and module == @module')['submodules']
            if len(submodules) == 0 or submodules.isna().all():
                submodules = _get_or_default_from_config(config['DATASETS'], config['defaults'], dataset, module)
            else:
                submodules = submodules.to_list()[0]
            if len(submodules) == 0:
                raise ValueError(
                    f'No submodules defined for {dataset}, {module}.'
                    ' Check that they are either defined either in the config or modules.tsv'
                )
            return submodules

        config['DATASETS'][dataset]['integration'] = _get(params, dataset, 'integration')
        config['DATASETS'][dataset]['metrics'] = _get(params, dataset, 'metrics')

    config = set_defaults(config, ['integration', 'metrics'])
    config['dataset_meta'] = str(Path(config['dataset_meta']).resolve())

    return config


def config_for_module(config, module):
    if 'DATASETS' not in config:
        config['DATASETS'] = {}
    datasets = config['DATASETS']
    for dataset in datasets.keys():
        # for module in datasets[dataset]['input'].keys():
        input_slot = datasets[dataset]['input']
        if input_slot is None or module not in input_slot:
            continue
        input_file = input_slot[module]
        if input_file in config['output_map']:
            file_name = config['output_map'][input_file].format(dataset=dataset)
            config['DATASETS'][dataset]['input'][module] = file_name
    return config


def update_input_files_per_dataset(
    dataset: str,
    module_name: str,
    config: dict,
    first_module: str = None,
    config_class_map: dict[str: ModuleConfig] = None,
    config_kwargs: dict[str: dict] = None,
):
    if config_kwargs is None:
        config_kwargs = {}
    if module_name == first_module:
        raise ValueError(f'Circle detected: first module {first_module} cannot be an input module')
    if first_module is None:
        first_module = module_name
    if config_class_map is None:
        config_class_map = {}
    
    file_map = InputFiles.parse(config['DATASETS'][dataset]['input'][module_name])
    input_files = {}
    for file_name, file_path in file_map.items():
        if '/' in file_path:
            # assert Path(file_path).exists(), f'Missing input file "{file_path}"'
            input_files |= {file_name: file_path}
            continue
        
        # get output files for input module
        input_module = file_path
        config = update_input_files_per_dataset(
            dataset=dataset,
            module_name=input_module,
            config=config,
            first_module=first_module,
            config_class_map=config_class_map,
        )
        
        ModuleConfigClass = config_class_map.get(input_module, ModuleConfig)
        input_cfg = ModuleConfigClass(
            module_name=input_module,
            config=config,
            datasets=[dataset],
            **config_kwargs.get(input_module, {})
        )
        input_files |= InputFiles.parse(
            input_cfg.get_output_files(
                subset_dict={'dataset': dataset},
                as_dict=True
            )
        )
    
    config['DATASETS'][dataset]['input'][module_name] = input_files
    return config
