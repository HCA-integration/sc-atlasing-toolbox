import yaml
import glob
import scanpy as sc


with open('test/config.yaml', 'r') as file:
    config = yaml.safe_load(file)


files = glob.glob('test/out/subset/*/*.h5ad')
print(files)

for file in files:
    adata = sc.read(file)
    assert 'subset' in adata.uns

    n_cells = adata.n_obs
    n_cells_threshold = config['DATASETS']['blood']['subset']['n_cells']

    print(file)
    print('\tthreshold: ', n_cells_threshold)
    print('\tadata.n_obs: ', n_cells)
    assert adata.n_obs < n_cells_threshold
