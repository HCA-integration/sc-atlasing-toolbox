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
    
    print(f'read {file}...')
    with zarr.open(file) as z:
        obs = read_elem(z["obs"])
    
    # new column mapping
    new_label_cfg = config['DATASETS'][dataset]['relabel'].get('new_columns')
    mapping_df = pd.read_table(new_label_cfg['file'])

    for label in new_label_cfg['order']:
        assert label in obs.columns
        assert label in mapping_df.columns
        print(obs[label].dtype)
        print(obs[label].value_counts())
    
    # existing column merge
    merge_cfg = config['DATASETS'][dataset]['relabel'].get('merge_columns')
    if merge_cfg is None:
        continue
    sep = merge_cfg.get('order')
    merge_df = pd.read_table(merge_cfg['file'])
    # TODO: check
