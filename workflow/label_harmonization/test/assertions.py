import anndata
import glob

files = glob.glob('test/out/label_harmonisation/*/*.zarr')
if len(files) == 0:
    print('No files found to test')

for file in files:
    print(f'read {file}...')
    adata = anndata.read_zarr(file)
    print(adata.obs['cell_type'].value_counts())
    print(adata.obs['harmonized_label'].value_counts())
    print(adata.obs['lineage'].value_counts())
