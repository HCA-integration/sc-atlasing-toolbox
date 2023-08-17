from anndata.experimental import read_elem
import zarr
import glob
import yaml
import pandas as pd


files = glob.glob('test/out/relabel/*.zarr')
if not files:
    print('No files found to test')

config = yaml.safe_load(open('test/config.yaml', 'r'))

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    mapping_config = config['DATASETS'][dataset]['relabel']['mapping']

    print(f"read {mapping_config['file']}...")
    mapping_df = pd.read_table(mapping_config['file'])

    print(f'read {file}...')
    z = zarr.open(file)
    obs = read_elem(z["obs"])

    for label in mapping_config['order']:
        assert label in obs.columns
        assert label in mapping_df.columns
        print(obs[label].dtype)
        print(obs[label].value_counts())
    
