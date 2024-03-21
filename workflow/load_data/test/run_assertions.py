from pprint import pprint
import glob
from anndata.experimental import read_elem
import zarr

from load_data_utils import SCHEMAS

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
    for col in SCHEMAS['CELLxGENE_OBS'] + SCHEMAS['EXTRA_COLUMNS']:
        assert col in obs.columns, f'"{col}" not in obs for "{file}"'
        assert not obs[col].isna().all(), f'"{col}" is all NA in obs for "{file}"'

    for col in SCHEMAS['CELLxGENE_VARS']:
        assert col in var.columns, f'"{col}" not in var for "{file}"'
        assert not var[col].isna().all(), f'"{col}" is all NA in var for "{file}"'

    for col in ['organ', 'dataset']:
        assert col in uns, f'"{col}" not in uns for "{file}"'

