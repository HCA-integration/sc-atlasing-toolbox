from utils.wildcards import wildcards_to_str
from utils.misc import ifelse


rule embedding:
    input:
        anndata='{filename}.h5ad'
    output:
        plot='{filename}_embedding.png'
    params:
        color='bulk_labels',
        basis='X_pca',
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/embedding.py'


rule umap:
    input:
        anndata='{filename}.h5ad'
    output:
        plot='{filename}_umap.png',
        coordinates='{filename}_coordinates.npy',
    params:
        color='bulk_labels',
        use_rep='X_pca',
    conda:
        ifelse(
            'os' not in config.keys() or config['os'] == 'm1',
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    script:
        '../scripts/umap.py'


rule barplot:
    input:
        tsv='test/data/integration.benchmark.tsv'
    output:
        png='{out_dir}/barplot_{metric}.png',
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        hue=None,
        facet_row='dataset',
        facet_col=None,
        title='Barplot',
        description=wildcards_to_str,
        dodge=True,
        xlim=(-.01, None),
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/barplot.py'


rule swarmplot:
    input:
        tsv='test/data/metrics.tsv'
    output:
        png='{out_dir}/swarmplot_{metric}.png',
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        title='Swarmplot',
        description=wildcards_to_str,
        ylim=(-.05, 1.05),
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/swarmplot.py'
