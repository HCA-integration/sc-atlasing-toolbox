from utils.environments import get_env


rule run_method:
    input:
        h5ad='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}/{method}/adata.zarr'),
        model=touch(directory('{dataset}/{method}/model'))
    params:
        norm_counts='normcounts',
        raw_counts='X',
        output_type='embed',
        hyperparams=None,
        env='scib',
    conda:
        lambda wildcards, params: get_env(config, params.env)
    # shadow: 'minimal'
    script:
        '../scripts/methods/{wildcards.method}.py'
