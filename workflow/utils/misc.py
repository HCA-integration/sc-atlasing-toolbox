import pandas as pd
import typing


def all_but(_list, is_not):
    return [x for x in _list if x != is_not]


def unique_dataframe(df):
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
