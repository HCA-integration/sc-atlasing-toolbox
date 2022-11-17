import pandas as pd
import scanpy as sc
import glob

from utils import CELLxGENE_OBS, CELLxGENE_VARS

datasets_df = pd.read_table('test/datasets.tsv')

# single dataset
single_outputs = glob.glob('test/out/*/processed/*.h5ad')
print(single_outputs)

for file in single_outputs:
    print(f'Check {file}...')
    adata = sc.read(file)

    assert 'dataset' in adata.uns
    assert 'dataset' in adata.obs
    assert 'organ' in adata.obs
    assert 'meta' in adata.uns
    for key in datasets_df.columns:
        assert key in adata.uns['meta']
        assert key in adata.obs.columns

# merged datasets

merged_outputs = glob.glob('test/out/*/merged/*.h5ad')

for file in merged_outputs:
    print(f'Check {file}...')
    adata = sc.read(file)
    print('Check that CELLxGENES mandatory columns present')
    extra_columns = ['organ', 'donor', 'sample', 'cell_annotation', 'reference', 'study', 'dataset', 'dataset_id']
    try:
        for col in CELLxGENE_OBS + extra_columns:
            assert col in adata.obs.columns

        for col in CELLxGENE_VARS:
            assert col in adata.var.columns

        for col in ['organ', 'dataset']:
            assert col in adata.uns
    except AssertionError:
        print(f'column {col} not in adata')
    finally:
        print(adata)
