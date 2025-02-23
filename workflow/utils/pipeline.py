"""
Utils that are specific to the overall pipeline but not the modules
"""

from pathlib import Path

from .config import _get_or_default_from_config, set_defaults
from .ModuleConfig import ModuleConfig
from .InputFiles import InputFiles


# deprecated
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


# deprecated
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
    if module_name == first_module:
        raise ValueError(f'Circle detected: first module {first_module} cannot be an input module')
    if first_module is None:
        first_module = module_name
    if config_class_map is None:
        config_class_map = {}
    if config_kwargs is None:
        config_kwargs = {}
    
    file_map = InputFiles.parse(config['DATASETS'][dataset]['input'][module_name])
    input_files = {}
    for file_name, file_path in file_map.items():
        if ',' in file_path and '.zarr' not in file_path and '.h5ad' not in file_path:
            # multiple files
            file_paths = file_path.split(',')
        else:
            file_paths = [file_path]
        
        for file_path in file_paths:
            file_path = file_path.strip()
            if '/' in file_path:
                # assert Path(file_path).exists(), f'Missing input file "{file_path}"'
                input_files |= {file_name: file_path}
                continue
            
            # get output files for input module
            input_module = file_path  # rename for easier readability
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
                warn=False,
                **config_kwargs.get(input_module, {})
            )
            output_files = input_cfg.get_output_files(
                subset_dict={'dataset': dataset},
                as_dict=True
            )
            input_files |= InputFiles.parse(output_files)
    
    config['DATASETS'][dataset]['input'][module_name] = input_files
    return config


def update_file_for_module_param(
    dataset: str,
    module_name: str,
    config: dict,
    key: str,
    subset_dict: dict = None,
    config_class_map: dict[str: ModuleConfig] = None,
    config_kwargs: dict[str: dict] = None,
):
    """
    :param dataset: dataset name
    :param module_name: module name
    :param config: config dict
    :param key: key in module that should specify a file name
    """
    if config_class_map is None:
        config_class_map = {}
    if config_kwargs is None:
        config_kwargs = {}
    
    dataset_config = config['DATASETS'].get(dataset, {})
    file_pattern = dataset_config.get(module_name, {}).get(key)
    if file_pattern in dataset_config['input']:
        input_module = file_pattern
        ModuleConfigClass = config_class_map.get(input_module, ModuleConfig)
        input_cfg = ModuleConfigClass(
            module_name=input_module,
            config=config,
            **config_kwargs.get(input_module, {})
        )
        if subset_dict is None:
            subset_dict = {}
        subset_dict |= {'dataset': dataset}
        dataset_config[module_name][key] = input_cfg.get_output_files(subset_dict=subset_dict)[0]
    return config