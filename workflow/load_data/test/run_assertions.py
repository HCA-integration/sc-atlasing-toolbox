import pandas as pd
import scanpy as sc
import glob

from utils import CELLxGENE_COLUMNS


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

# merged datasets

merged_outputs = glob.glob('test/out/*/merged/*.h5ad')

for file in merged_outputs:
    print(f'Check {file}...')
    adata = sc.read(file)

    print('Check that CELLxGENES mandatory columns present')
    for col in CELLxGENE_COLUMNS:
        assert col in adata.obs.columns

    assert 'organ' in adata.uns
    assert 'dataset' in adata.uns
    assert 'organ' in adata.obs
    assert 'dataset' in adata.obs
    assert 'donor' in adata.obs
    assert 'sample' in adata.obs
    print(adata)