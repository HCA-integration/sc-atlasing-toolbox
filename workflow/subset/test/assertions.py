import scanpy as sc
import glob

files = glob.glob('test/out/subset/*/*.h5ad')
print(files)

for file in files:
    adata = sc.read(file)
    assert 'subset' in adata.uns

