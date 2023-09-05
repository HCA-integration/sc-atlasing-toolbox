import yaml
import glob
from anndata.experimental import read_elem
import zarr


with open('test/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

files = glob.glob('test/out/subset/*/*.zarr')
# print(files)

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    n_cells_threshold = config['DATASETS'][dataset]['subset']['n_cells']

    print(f'read {file}...')
    z = zarr.open(file)
    obs = read_elem(z["obs"])
    uns = read_elem(z["uns"])

    print(dataset)
    assert 'subset' in uns
    n_cells = obs.shape[0]
    print('\tthreshold: ', n_cells_threshold)
    print('\tadata.n_obs: ', n_cells)
    assert n_cells < n_cells_threshold
