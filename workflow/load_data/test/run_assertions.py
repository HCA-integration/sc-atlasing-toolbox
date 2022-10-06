import scanpy as sc

file = 'test/out/download/SciTranslMed2021_covid19.h5ad'
adata = sc.read(file)

# TODO assertions