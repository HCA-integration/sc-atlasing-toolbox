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
        tsv=expand('{resolution}.tsv', resolution=[0.5, 1.0]),
    output:
        tsv='all_resolutions.tsv'
    run:
        from utils.misc import merge

        dfs = [pd.read_table(file, index_col=0) for file in input.tsv]
        cluster_df = merge(
            dfs,
            left_index=True,
            right_index=True,
        )
        print(cluster_df)
        cluster_df.to_csv(output.tsv, sep='\t')


rule plot_umap:
    input:
        zarr='dataset.h5ad',
        clusters='all_resolutions.tsv'
    output:
        png='clusters.png',
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot_umap.py'