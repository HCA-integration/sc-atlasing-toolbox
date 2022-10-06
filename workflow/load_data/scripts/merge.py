import scanpy as sc

files = snakemake.input
out_file = snakemake.output.h5ad
print(files)

adatas = [sc.read(file) for file in files]
print(adatas)

adata = sc.concat(adatas)
adata.write(out_file)