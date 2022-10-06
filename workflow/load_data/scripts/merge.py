import scanpy as sc

organ = snakemake.wildcards.organ
files = snakemake.input
out_file = snakemake.output.h5ad

adatas = [sc.read(file) for file in files]
print(adatas)

print(f'Concatenate...')
adata = sc.concat(adatas)

print(f'Add metadata')
adata.uns['dataset'] = organ
adata.uns['organ'] = organ
adata.obs['organ'] = organ

adata.write(out_file)