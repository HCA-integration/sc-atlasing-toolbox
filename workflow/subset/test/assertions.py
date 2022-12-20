import yaml
import anndata
import glob


with open('test/config.yml', 'r') as file:
    config = yaml.safe_load(file)


files = glob.glob('test/out/subset/*/*.zarr')
print(files)

for file in files:
    adata = anndata.read_zarr(file)
    assert 'subset' in adata.uns
    print(adata.n_obs)
    assert adata.n_obs < config['subset']['blood']['n_cells']
