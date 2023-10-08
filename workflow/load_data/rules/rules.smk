from utils.environments import get_env


rule merge:
    input:
        [
            'dataset1.h5ad',
            '{dataset2.h5ad'
        ]
    output:
        zarr=directory('merged.zarr'),
    params:
        dataset='dataset',
        merge_strategy='inner'
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    script:
        '../scripts/merge.py'


rule filter:
    """
    Apply QC filters
    """
    input:
        zarr='merged.zarr'
    output:
        zarr=directory('filtered.zarr'),
        removed=directory('removed.zarr'),
    params:
        filter={
            'cells_per_sample': {'min': 50, 'max': 10000},
            'mito_pct': 30,
            'remove_by_colum': {'obs_column': ['entry_to_exclude_1', 'entry_to_exclude_2']}
        }
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    script:
        '../scripts/filter.py'
