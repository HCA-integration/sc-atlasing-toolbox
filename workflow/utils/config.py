import warnings
import pandas as pd

from .misc import expand_dict


def set_defaults(config, modules=None, warn=True):
    if 'defaults' not in config:
        config['defaults'] = {}
    if 'datasets' not in config['defaults']:
        config['defaults']['datasets'] = list(config['DATASETS'].keys())

    if modules is None:
        modules = ['integration', 'metrics']
    elif isinstance(modules, str):
        modules = [modules]

    for module in modules:
        for dataset in config['DATASETS'].keys():
            entry = _get_or_default_from_config(
                config['DATASETS'],
                config['defaults'],
                dataset,
                module,
                warn=warn,
            )
            # for TSV input make sure integration methods have the proper types
            if module == 'integration' and isinstance(entry, list):
                # get parameters from config
                entry = {k: config['defaults'][module][k] for k in entry}
            config['DATASETS'][dataset][module] = entry
    return config


def get_wildcards_from_config(
        config,
        config_params,
        wildcard_names,
        explode_by=None,
        config_keys=None,
):
    """

    :param config: Part of the Snakemake config dictionary. The config_params must be contained per entry.
    :param config_params: list of parameters for each config entry
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param explode_by: column to explode by, expecting list entry for that column
    :param config_keys: list of entries to subset the config by., otherwise use all keys
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """

    if config_keys is None:
        config_keys = config.keys()
    records = [
        (key, *[_get_or_default_from_config(config, {}, key, w) for w in config_params])
        for key in config_keys
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df.reset_index(drop=True)


def _get_or_default_from_config(config, defaults, key, value, warn=True):
    """
    Get entry from config or return defaults if not present

    :param config: part of the config with multiple entries of the same structure
    :param defaults: config defaults for keys with missing value
    :param key: top-level key of config dictionary
    :param value: key of
    :return:
    """
    if key not in config.keys():
        print('config:', config)
        raise KeyError(f'Key "{key}" not found in config')

    if value in config[key]:
        return config[key][value]
    try:
        assert value in defaults.keys()
    except AssertionError:
        if warn:
            warnings.warn(f'No default defined for "{value}" for "{key}", returning None')
        return None
    return defaults[value]


def get_hyperparams(config, module='integration'):
    """
    Get hyperparameters specific to each method of a module for all datasets

    :param config: config containing dataset specific information
    :param module: name of module, key must be present for each dataset entry
    :return: DataFrame with hyperparameters
    """
    records = []
    for dataset, dataset_dict in config['DATASETS'].items():
        for method, hyperparams_dict in dataset_dict[module].items():
            if isinstance(hyperparams_dict, dict):
                for rec in expand_dict(hyperparams_dict):
                    records.append((dataset, method, *rec))
            else:
                records.append((dataset, method, str(hyperparams_dict), hyperparams_dict))
    return pd.DataFrame(records, columns=['dataset', 'method', 'hyperparams', 'hyperparams_dict'])


def get_resource(config, resource_key, profile='cpu'):
    """
    Retrieve resource information from config['resources']
    
    :param config: config passed from Snakemake
    :param profile: resource profile, key under config['resources']
    :param resource_key: resource key, key under config['resources'][profile]
    """
    if 'resources' not in config or not profile:
        return ''
    resources = config['resources']
    try:
        res = resources[profile][resource_key]
    except KeyError:
        print(
            f'WARNING: Invalid profile "{profile}" or resource key "{resource_key}". '
            'Please check that your config contains the correct entries under config["resources"]'
        )
        res = ''
    return res


def get_datasets_for_module(config, module):
    """
    Collect dataset names e.g. for wildcard expansion from config["DATASETS"] for a given module
    If the "DATASETS" key is not available in the config, warn and return an empty list
    A dataset is valid if it contains an input file for the given module.

    :param config: config dictionary passed from Snakemake
    :param module: name of module to collect valid datasets for
    :return: list of valid dataset names from config["DATASETS"]
    """
    if 'DATASETS' not in config:
        warnings.warn('No datasets specified in config, cannot collect any datasets')
        return []
    return [
        dataset for dataset in config['DATASETS'].keys()
        if 'input' in config['DATASETS'][dataset]
           and config['DATASETS'][dataset]['input'] is not None
           and module in config['DATASETS'][dataset]['input']
    ]