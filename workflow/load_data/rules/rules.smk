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
