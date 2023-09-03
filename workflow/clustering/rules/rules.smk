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