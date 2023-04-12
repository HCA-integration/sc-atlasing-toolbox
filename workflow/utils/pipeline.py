"""
Utils that are specific to the overall pipeline but not the modules
"""

from pathlib import Path

from .config import _get_or_default_from_config, set_defaults


def update_module_configs(config, params):
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
