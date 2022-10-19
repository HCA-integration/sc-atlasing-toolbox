import numpy as np
import pandas as pd


def get_wildcards_from_config(config, config_params, wildcard_names, explode_by, config_keys=None):
    """

    :param config: part of the config with entries of same format
    :param config_params: list of parameters for each config entry
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param config_keys: list of entries to subset the config by., otherwise use all keys
    :param explode_by: column to explode by, expecting list entry for that column
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """
    if config_keys is None:
        config_keys = config.keys()
    records = [
        (key, *[config[key][w] for w in config_params])
        for key in config_keys
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df.reset_index(drop=True)


def get_params(wildcards, parameters_df, column, wildcards_keys=None):
    """

    :param wildcards: wildcards passed on from Snakemake or dictionary
    :param parameters_df: dataframe with parameters for wildcards, column names must match wildcard names
    :param column: column or columns from parameters_df
    :param wildcards_keys: list of wildcards to constrain the query to
    :return: single parameter or list of parameters as specified by column
    """
    assert column in parameters_df.columns

    if wildcards_keys is None:
        wildcards_keys = [key for key in wildcards.keys() if key in parameters_df.columns]

    assert np.all([key in wildcards.keys() for key in wildcards.keys()])

    query = ' and '.join([f'{w} == @wildcards.{w}' for w in wildcards_keys])
    params_sub = parameters_df.query(query)
    try:
        assert params_sub.shape[0] == 1
    except AssertionError:
        raise ValueError(f'Incorrect number of parameters for wildcards: {wildcards}\n{parameters_df}')
    param = params_sub[column].tolist()
    if len(param) == 1:
        return param[0]
    return param


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
