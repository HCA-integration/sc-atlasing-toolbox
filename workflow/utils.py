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
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df


def get_params(wildcards, parameters_df, column):
    """

    :param wildcards: wildcards passed on from Snakemake
    :param parameters_df: dataframe with parameters for wildcards, column names must match wildcard names
    :param column: column or columns from parameters_df
    :return: single parameter or list of parameters as specified by column
    """
    assert column in parameters_df.columns
    query = ' and '.join([f'{w} == @wildcards.{w}' for w in wildcards.keys()])
    params_sub = parameters_df.query(query)
    param = params_sub[column].unique().tolist()
    if len(param) == 1:
        return param[0]
    return param