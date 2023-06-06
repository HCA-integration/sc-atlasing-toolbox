import numpy as np
import pandas as pd
import typing


def remove_outliers(adata, extrema='max', factor=10):
    umap = adata.obsm['X_umap']
    if extrema == 'max':
        abs_values = np.abs(umap.max(axis=1))
    elif extrema == 'min':
        abs_values = np.abs(umap.min(axis=1))
    outlier_mask = abs_values < factor * abs_values.mean()
    return adata[outlier_mask]


def all_but(_list, is_not):
    return [x for x in _list if x != is_not]


def unique_dataframe(df):
    if df.empty:
        return df
    hashable_columns = [col for col in df.columns if isinstance(df[col].iloc[0], typing.Hashable)]
    duplicated = df[hashable_columns].duplicated()
    return df[~duplicated]


def expand_dict(_dict):
    """
    Create a cross-product on a dictionary with literals and lists
    :param _dict: dictionary with lists and literals as values
    :return: list of dictionaries with literals as values
    """
    df = pd.DataFrame({k: [v] if isinstance(v, list) else [[v]] for k, v in _dict.items()})
    for col in df.columns:
        df = df.explode(col)
    dict_list = df.apply(lambda row: {col: x for col, x in zip(df.columns, row)}, axis=1)

    def remove_chars(s, chars="{} ',"):
        for c in chars:
            s = s.replace(c, '')
        return s

    wildcards = df.apply(lambda row: '-'.join([remove_chars(f'{col}:{x}') for col, x in zip(df.columns, row)]), axis=1)
    return zip(wildcards, dict_list)


def unlist_dict(dictionary):
    return {
        k: v[0] if isinstance(v, list) and len(v) == 1 else v
        for k, v in dictionary.items()
    }

def unpack_dict_in_df(df, column):
    """
    Given a column in a pandas dataframe containing dictionies, extract these to top level
    
    :param df: pandas dataframe
    :param column: column name containing dictionaries 
    """
    return df.drop(columns=column).assign(**df[column].dropna().apply(pd.Series, dtype=object))


def ifelse(statement, _if, _else):
    if statement:
        return _if
    else:
        return _else
