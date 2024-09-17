import pandas as pd


ALL_SLOTS = [
    'X',
    'layers',
    'raw',
    'obs',
    'obsm',
    'obsp',
    'var',
    'varm',
    'varp',
    'uns',
]

def get_same_columns(adatas):
    _ad1 = next(iter(adatas.values()))
    # collect obs columns that exist in all files
    same_obs_columns = set(
        _ad1.obs.columns.tolist()
    ).intersection(
        *[_ad.obs.columns.tolist() for _ad in adatas.values()]
    )
    # determine which columns do not have the same values across all files
    obs_to_remove = []
    for col in same_obs_columns:
        same_across = all(
            _ad1.obs[col].sort_values().equals(_ad.obs[col].sort_values())
            for _ad in adatas.values()
        )
        if not same_across:
            obs_to_remove.append(col)
    # remove columns that are not the same across all files
    return [
        col for col in same_obs_columns
        if col not in obs_to_remove
    ]


def merge_df(
    df_current,
    file_id,
    df_previous,
    same_columns,
    sep='_',
    obs_index_col=None,
):
    if obs_index_col is not None:
        assert obs_index_col in df_current.columns, f'Index column "{obs_index_col}" not found for {file_id}\n{df_current}'
        df_current = df_current.set_index(obs_index_col)
    if df_previous is None:
        return df_current.rename(
            columns={
                col: f'{col}{sep}{file_id}'
                for col in df_current.columns.tolist()
                if col not in same_columns
            }
        )
    
    idx_diff = df_current.index.difference(df_previous.index).sort_values()
    n_diff = len(idx_diff)
    assert n_diff == 0, \
        f'Index must be the same\n {n_diff} differing indices: {idx_diff}\n' \
        f'current index: {df_current.index.sort_values()} \nprevious index: {df_previous.index.sort_values()}'
    unique_columns = [col for col in df_current.columns if col not in same_columns]
    df_current = df_current[unique_columns].rename(
        columns={col: f'{col}{sep}{file_id}' for col in unique_columns}
    )
    df = pd.concat([df_previous, df_current], axis=1)
    return df
