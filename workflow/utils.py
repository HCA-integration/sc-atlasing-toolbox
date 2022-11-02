import numpy as np
import pandas as pd
from snakemake.io import expand


def set_defaults(config):
    if 'defaults' not in config:
        config['defaults'] = {}
    if 'datasets' not in config['defaults']:
        config['defaults']['datasets'] = list(config['DATASETS'].keys())
    return config


def get_wildcards_from_config(
        config,
        config_params,
        wildcard_names,
        explode_by,
        config_keys=None,
        slot='DATASETS'
):
    """

    :param config: Snakemake config dictionary
    :param config_params: list of parameters for each config entry
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param explode_by: column to explode by, expecting list entry for that column
    :param config_keys: list of entries to subset the config by., otherwise use all keys
    :param slot: slot in config to get wildcards from. The config_params must be contained under that slot.
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """
    global_config = config
    config = config[slot]

    if config_keys is None:
        config_keys = config.keys()
    records = [
        (key, *[_get_or_default_from_config(config, global_config['defaults'], key, w) for w in config_params])
        for key in config_keys
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df.reset_index(drop=True)


def _get_or_default_from_config(config, defaults, key, value):
    """
    Get entry from config or return defaults if not present

    :param config: part of the config with multiple entries of the same structure
    :param defaults: config defaults for keys with missing value
    :param key: top-level key of config dictionary
    :param value: key of
    :return:
    """
    if value in config[key]:
        return config[key][value]
    try:
        assert value in defaults
    except AssertionError:
        raise AssertionError(f'No default defined for "{value}"')
    return defaults[value]


def get_wildcards(wildcards_df, columns, wildcards=None):
    """
    Get wildcards from DataFrame

    :param wildcards_df: DataFrame with wildcard names in columns
    :param columns: wildcard keys that are present in the dataframe column
    :param wildcards: wildcards passed from Snakemake to subset to
    :return: subset of the wildcards_df by wildcard match and columns
    """
    if wildcards:
        query = ' and '.join([f'{w} == @wildcards.{w}' for w in wildcards.keys()])
        wildcards_df = wildcards_df.query(query)
    return wildcards_df[columns].drop_duplicates().to_dict('list')


def expand_per(target, wildcards_df, wildcards, columns):
    """
    Expand a target by specific wildcards, considering wildcards to

    :param target: file name or Snakemake target object with wildcards
    :param wildcards_df: DataFrame with wildcard names in columns
    :param wildcards: wildcards passed from Snakemake to expand
    :param columns: wildcard keys that are present in the dataframe column
    :return: list of expanded targets
    """
    return expand(target, zip, **get_wildcards(wildcards_df, columns, wildcards), allow_missing=True)


def get_params(wildcards, parameters_df, column, wildcards_keys=None):
    """

    :param wildcards: wildcards passed on from Snakemake or dictionary
    :param parameters_df: dataframe with parameters for wildcards, column names must match wildcard names
    :param column: column or columns from parameters_df
    :param wildcards_keys: list of wildcards to return
    :return: single parameter or list of parameters as specified by column
    """
    try:
        assert column in parameters_df.columns

        if wildcards is not None:
            if wildcards_keys is None:
                wildcards_keys = [key for key in wildcards.keys() if key in parameters_df.columns]
            assert np.all([key in wildcards.keys() for key in wildcards.keys()])

            query = ' and '.join([f'{w} == @wildcards.{w}' for w in wildcards.keys()])
            params_sub = parameters_df.query(query)
        else:
            params_sub = parameters_df

        try:
            assert params_sub.shape[0] == 1
        except AssertionError:
            raise ValueError(
                f'Incorrect number of parameters\n{params_sub}'
            )
        param = params_sub[column].tolist()
    except Exception as e:
        raise Exception(
            f'Error for wildcards={wildcards}, column={column}, wildcards_keys={wildcards_keys}'
            f'\n{parameters_df}\nError message: {e}'
        )
    if len(param) == 1:
        return param[0]
    return param[wildcards_keys]


def get_resource(config, profile, resource_key):
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
    except:
        raise KeyError(
            f'Invalid profile "{profile}" or resource key "{resource_key}". '
            'Please check that your config contains the correct entries under config["resources"]'
        )
    return res
