"""
Preprocessing steps
Only the unique outputs per step are saved for storage efficiency. For assembled zarr files, see `assemble.smk`.
"""
from utils.misc import ifelse


rule normalize:
    input: '{dataset}.h5ad'
    output:
        zarr=directory('{dataset}_normalized.zarr'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/normalize.py'


rule highly_variable_genes:
    input:
        zarr='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_highly_variable_genes.zarr')
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/highly_variable_genes.py'


rule pca:
    input:
        zarr='{dataset}.h5ad',
        counts='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_pca.zarr')
    conda: '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/pca.py'


rule neighbors:
    input:
        zarr='{dataset}.h5ad'
    output:
        zarr=directory('{dataset}_neighbors.zarr')
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/neighbors.py'


rule umap:
    input:
        zarr='{dataset}.h5ad',
        rep='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_umap.zarr')
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/umap.py'


### Assemble ###

rule assemble:
    input:
        unpack(
            dict(
                counts='{dataset}.h5ad',
                normalize='{dataset}_normalized.zarr',
                highly_variable_genes='{dataset}_highly_variable_genes.zarr',
                pca='{dataset}_pca.zarr',
                neighbors='{dataset}_neighbors.zarr',
                umap='{dataset}_umap.zarr',
            )
        )
    output:
        zarr=directory('{dataset}_preprocessed.zarr')
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/assemble.py'


### Plots ###

rule plot_embedding:
    input:
        anndata='{dataset}.h5ad',
    output:
        plot='plot_{dataset}_pca.png',
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/plot_embedding.py'


rule plot_umap:
    input:
        anndata='{dataset}.h5ad'
    output:
        plot='plot_{dataset}.png',
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    script:
        '../scripts/plot_umap.py'
