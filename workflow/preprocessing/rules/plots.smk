from utils.wildcards import wildcards_to_str
from utils.misc import ifelse


rule plot_embedding:
    input:
        anndata=rules.pca.output.zarr,
    output:
        plot=image_dir / 'pca' / '{dataset}.png',
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'colors']),
        basis='X_pca',
        ncols=1,
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/plot_embedding.py'


rule plot_umap:
    input:
        anndata=rules.umap.output.zarr
    output:
        plot=image_dir / 'umap' / '{dataset}.png',
    params:
        color=lambda w: get_for_dataset(config, w.dataset, [module_name, 'colors']),
        ncols=1,
        outlier_factor=3
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    script:
        '../scripts/plot_umap.py'
