from utils.environments import get_env


rule filter:
    input:
        zarr='dataset.zarr'
    output:
        zarr=directory('dataset_filtered.zarr'),
    conda:
        get_env(config, 'scanpy', env_dir='envs')
    script:
        '../scripts/filter.py'
