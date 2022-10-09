import pandas as pd


def get_wildcards_from_config(config, config_params, wildcard_names, explode_by):
    """

    :param config: part of the config with entries of same format
    :param config_params: list of parameters for each config entry
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param explode_by: column to explode by, expecting list entry for that column
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """
    records = [
        (key, *[config[key][w] for w in config_params])
        for key in config.keys()
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        df = df.explode(explode_by)
    return df


def get_params(wildcards, parameters, column):
    assert column in parameters.columns
    query = ' and '.join([f'{w} == @wildcards.{w}' for w in wildcards.keys()])
    params_sub = parameters.query(query)
    param = params_sub[column].unique().tolist()
    if len(param) == 1:
        return param[0]
    return param