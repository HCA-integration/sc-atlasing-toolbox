import scanpy as sc

in_file = snakemake.input.h5ad
out_file = snakemake.output.h5ad
wildcards = snakemake.wildcards
meta = snakemake.params.meta

adata = sc.read(in_file)
print(adata)

adata.uns['dataset'] = wildcards.dataset
adata.obs['dataset'] = wildcards.dataset
adata.obs['organ'] = meta['organ']
adata.uns['meta'] = meta

# TODO plot count distribution -> save to file

adata.write(out_file)