import pandas as pd
from utils.environments import get_env


rule cluster:
    input:
        zarr='dataset.h5ad'
    output:
        tsv='{resolution}.tsv',
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    script:
        '../scripts/clustering.py'


rule merge:
    input:
        zarr='dataset.h5ad',
        tsv=expand('{resolution}.tsv', resolution=[0.5, 1.0]),
    output:
        zarr=directory('all_resolutions.zarr')
    localrule: True
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/merge.py'