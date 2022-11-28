import scanpy as sc
import glob

# celltypist
celltypist = glob.glob('test/out/label_transfer/celltypist/*/*.h5ad')
print(celltypist)

for file in celltypist:
    adata = sc.read(file)
    celltypist_columns = [x for x in adata.obs.columns if x.startswith('celltypist')]
    print(celltypist_columns)
    try:
        assert len(celltypist_columns) == 4
    except AssertionError:
        print(adata)
