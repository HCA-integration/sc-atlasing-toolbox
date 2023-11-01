import numpy as np
import pandas as pd
import typing
import hashlib


def remove_outliers(adata, extrema='max', factor=10, rep='X_umap'):
    umap = adata.obsm[rep]
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
    # hashable_columns = [
    #     col for col in df.columns
    #     if all(isinstance(df[col].iloc[i], typing.Hashable) for i in range(df.shape[0]))
    # ]
    # duplicated = df[hashable_columns].duplicated()
    duplicated = df.astype(str).duplicated()
    return df[~duplicated].reset_index(drop=True)


def expand_dict(_dict):
    """
    Create a cross-product on a dictionary with literals and lists
    :param _dict: dictionary with lists and literals as values
    :return: list of dictionaries with literals as values
    """
    df = pd.DataFrame({k: [v] if isinstance(v, list) else [[v]] for k, v in _dict.items()})
    for col in df.columns:
        df = df.explode(col)
    dict_list = df.apply(lambda row: dict(zip(df.columns, row)), axis=1)

    def remove_chars(s, chars="{} ',"):
        for c in chars:
            s = s.replace(c, '')
        return s

    wildcards = df.apply(
        lambda row: '-'.join(
            [
                remove_chars(f'{col[0]}:{x}') for col, x in zip(df.columns, row)
            ]
        ),
        axis=1
    )
    return zip(wildcards, dict_list)


def expand_dict_and_serialize(_dict):
    """
    Create a cross-product on a dictionary with literals and lists
    :param _dict: dictionary with lists and literals as values
    :return: list of dictionaries with literals as values
    """
    import jsonpickle
    import hashlib

    df = pd.DataFrame({k: [v] if isinstance(v, list) else [[v]] for k, v in _dict.items()})
    for col in df.columns:
        df = df.explode(col)
    dict_list = df.apply(lambda row: dict(zip(df.columns, row)), axis=1)

    wildcards = [
        hashlib.blake2b(jsonpickle.encode(d).encode('utf-8'), digest_size=5).hexdigest()
        for d in dict_list
    ]

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


def ensure_sparse(adata):
    from scipy.sparse import csr_matrix, issparse

    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)


def merge(dfs, **kwargs):
    """
    Merge list of dataframes
    :param dfs: list of dataframes
    :param kwargs: arguments passed to pd.merge
    :return: merged dataframe
    """
    from functools import reduce

    merged_df = reduce(
        lambda x, y: pd.merge(x, y, **kwargs),
        dfs
    )
    print(merged_df)
    return merged_df


def create_hash(string: str, digest_size: int = 5):
    string = string.encode('utf-8')
    return hashlib.blake2b(string, digest_size=digest_size).hexdigest()