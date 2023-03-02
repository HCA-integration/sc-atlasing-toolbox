import yaml
import anndata
import glob


with open('test/config.yaml', 'r') as file:
    config = yaml.safe_load(file)


files = glob.glob('test/out/subset/*/*.zarr')
print(files)

for file in files:
    adata = anndata.read_zarr(file)
    assert 'subset' in adata.uns

    n_cells = adata.n_obs
    n_cells_threshold = config['DATASETS']['blood']['subset']['n_cells']

    print('adata.n_obs: ', n_cells)
    print('n_cells: ', n_cells_threshold)
    assert adata.n_obs < n_cells_threshold
