from pprint import pprint
import glob
from anndata.experimental import read_elem
import zarr

from utils import SCHEMAS

# single dataset
single_outputs = glob.glob('test/out/*/harmonize_metadata/*.zarr')
print('harmonised metadata output:')
pprint(single_outputs)

for file in single_outputs:
    print(f'Check {file}...')
    z = zarr.open(file)
    obs = read_elem(z["obs"])
    uns = read_elem(z["uns"])

    assert 'dataset' in uns
    assert 'dataset' in obs
    assert 'organ' in obs
    assert 'meta' in uns
    try:
        for key in SCHEMAS['CELLxGENE_OBS'] + SCHEMAS['EXTRA_COLUMNS']:
            assert key in obs.columns
    except AssertionError:
        raise AssertionError(f'Single dataset: "{key}" not in "{file}"')

# merged datasets

merged_outputs = glob.glob('test/out/*/merged/study/*.zarr') + glob.glob('test/out/*/merged/organ/*.zarr')
print('merged output:')
pprint( merged_outputs)

for file in merged_outputs:
    print(f'Check {file}...')
    z = zarr.open(file)
    obs = read_elem(z["obs"])
    uns = read_elem(z["uns"])
    var = read_elem(z["var"])
    print('Check that CELLxGENES mandatory columns present')
    try:
        for col in SCHEMAS['CELLxGENE_OBS'] + SCHEMAS['EXTRA_COLUMNS']:
            assert col in obs.columns
            assert not obs[col].isna().all()

        for col in SCHEMAS['CELLxGENE_VARS']:
            assert col in var.columns
            assert not var[col].isna().all()

        for col in ['organ', 'dataset']:
            assert col in uns
    except AssertionError:
        raise AssertionError(f'Merged dataset: column "{col}" not in "{file}"')

