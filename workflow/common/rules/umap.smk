rule umap:
    input:
        anndata='{filename}.h5ad'
    output:
        plot='{filename}.png'
    params:
        color='bulk_labels',
        use_rep='X_pca',
    conda:
        '../envs/scanpy_rapids.yaml'
    script:
        '../scripts/umap.py'