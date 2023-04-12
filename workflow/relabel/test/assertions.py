import scanpy as sc
import glob

files = glob.glob('test/out/label_harmonization/**/*.h5ad')
if len(files) == 0:
    print('No files found to test')

for file in files:
    print(f'read {file}...')
    adata = sc.read(file)

    for label in ['cell_type', 'harmonized_label', 'lineage']:
        print(adata.obs[label].dtype)
        print(adata.obs[label].value_counts())
    
