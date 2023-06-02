import pandas as pd
import anndata
import glob

from utils import SCHEMAS

# single dataset
single_outputs = glob.glob('test/out/*/processed/*.zarr')
print(single_outputs)

for file in single_outputs:
    print(f'Check {file}...')
    adata = anndata.read_zarr(file)

    assert 'dataset' in adata.uns
    assert 'dataset' in adata.obs
    assert 'organ' in adata.obs
    assert 'meta' in adata.uns
    try:
        for key in SCHEMAS['CELLxGENE_OBS'] + SCHEMAS['EXTRA_COLUMNS']:
            assert key in adata.obs.columns
    except AssertionError:
        raise AssertionError(f'Single dataset: "{key}" not in "{file}"')

# merged datasets

merged_outputs = glob.glob('test/out/*/merged/*.zarr')

for file in merged_outputs:
    print(f'Check {file}...')
    adata = anndata.read_zarr(file)
    print('Check that CELLxGENES mandatory columns present')
    try:
        for col in SCHEMAS['CELLxGENE_OBS'] + SCHEMAS['EXTRA_COLUMNS']:
            assert col in adata.obs.columns
            assert not adata.obs[col].isna().all()

        for col in SCHEMAS['CELLxGENE_VARS']:
            assert col in adata.var.columns
            assert not adata.var[col].isna().all()

        for col in ['organ', 'dataset']:
            assert col in adata.uns
    except AssertionError:
        raise AssertionError(f'Merged dataset: column "{col}" not in "{file}"')

