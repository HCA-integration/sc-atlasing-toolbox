import scanpy as sc

from utils import add_metadata, read_anndata, process

input_adata = snakemake.input.h5ad
output_adata = snakemake.output.h5ad
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = read_anndata(input_adata)
adata_raw.X = adata_raw.layers['normcounts'].copy()

# prepare output adata
adata = adata_raw
adata = process(adata=adata, adata_raw=adata_raw, output_type='full')
adata.obsm['X_emb'] = adata.obsm['X_pca']
sc.pp.neighbors(adata)
add_metadata(adata, wildcards, params)

adata.write(output_adata, compression='gzip')
