import anndata

input_files = snakemake.input
output_file = snakemake.output.h5ad

adatas = [anndata.read(file) for file in input_files]

adata = anndata.concat(adatas, merge='same', uns_merge='same')
print(adata)

adata.write(output_file)
