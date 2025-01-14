import pandas as pd
from utils.environments import get_env


rule cluster:
    input:
        zarr='dataset.h5ad'
    output:
        zarr='{resolution}.zarr',
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    script:
        '../scripts/clustering.py'


rule merge:
    input:
        zarr='dataset.h5ad',
        cluster_anno=expand('{resolution}.zarr', resolution=[0.5, 1.0]),
    output:
        zarr=directory('all_resolutions.zarr')
    localrule: True
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/merge.py'


rule plot_evaluation:
    input:
        zarr='dataset.h5ad',
    output:
        plots=directory('evaluation')
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot_enrichment.py'