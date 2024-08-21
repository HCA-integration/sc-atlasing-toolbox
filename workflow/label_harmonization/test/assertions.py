import scanpy as sc
import glob
from utils.io import read_anndata

files = glob.glob('test/out/label_harmonization/dataset~*/file_id~*/cellhint/*.h5ad')
if len(files) == 0:
    print('No files found to test')

for file in files:
    print(f'read {file}...')
    adata = read_anndata(file)

    for label in ['bulk_labels']: # ['cell_type', 'harmonized_label', 'lineage']:
        print(adata.obs[label].dtype)
        print(adata.obs[label].value_counts())
