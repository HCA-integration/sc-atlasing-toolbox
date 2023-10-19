"""
Helper functions for working with wildcards
"""

import traceback
import warnings

import numpy as np
import pandas as pd
from snakemake.io import expand

from .misc import unique_dataframe


def wildcards_to_str(wildcards):
    return ' '.join([f'{key}={value}' for key, value in wildcards.items()])


def get_wildcards(wildcards_df, columns=None, wildcards=None, drop_na=False):
    """
    Get wildcards from DataFrame

    :param wildcards_df: DataFrame with wildcard names in columns
    :param columns: wildcard keys that are present in the dataframe column
    :param wildcards: wildcards passed from Snakemake to subset to
    :param drop_na: whether to remove rows that contain NAs
    :return: subset of the wildcards_df by wildcard match and columns
    """
    if columns is None:
        columns = wildcards_df.columns
    if isinstance(columns, str):
        columns = [columns]
    wildcards_df = subset_by_wildcards(wildcards_df, wildcards)
    wildcards_df = unique_dataframe(wildcards_df[columns])
    if drop_na:
        wildcards_df = wildcards_df.dropna()
    return wildcards_df.to_dict('list')


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

        # quickfix: ignore hyperparameters
        if 'hyperparams' in wildcards_keys:
            wildcards_keys.remove('hyperparams')

        wildcards = {k: v for k, v in wildcards.items() if k in wildcards_keys}
        params_sub = subset_by_wildcards(parameters_df, wildcards)
        columns = list(set(wildcards_keys + [column]))
        params_sub = unique_dataframe(params_sub[columns])

        try:
            assert params_sub.shape[0] == 1
        except AssertionError:
            raise ValueError(
                f'More than 1 row after subsetting\n{params_sub}'
            )
        param = params_sub[column].tolist()
    except Exception as e:
        print(e)
        print(traceback.format_exc())
        raise Exception(
            f'Error for column={column}, wildcards={wildcards}, wildcards_keys={wildcards_keys}'
            f'\n{parameters_df}\nError message: {e}'
        )
    if len(param) == 1:
        return param[0]
    return param[wildcards_keys]


def subset_by_wildcards(df, wildcards):
    """
    Helper function to subset df by wildcards
    :param df: dataframe with wildcard names in column and wildcard values in rows
    :param wildcards: wildcards object from Snakemake
    :return: subset dataframe
    """
    if not wildcards:
        return df
    wildcards = {k: wildcards[k] for k in wildcards.keys() if k in df.columns}
    query = ' and '.join([f'{k} == "{v}"' for k, v in wildcards.items()])
    df = df.query(query)
    if df.shape[0] == 0:
        raise ValueError(f'no wildcard combination found in wildcards df {query}')
    return df
