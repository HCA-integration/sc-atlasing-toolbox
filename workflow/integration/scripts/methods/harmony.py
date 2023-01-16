from harmony import harmonize

from utils import add_metadata, read_anndata, process


input_adata = snakemake.input.h5ad
output_adata = snakemake.output.h5ad
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = read_anndata(input_adata)
adata_raw.X = adata_raw.layers['normcounts'].copy()

# run method
adata = adata_raw.copy()
adata.obsm["X_emb"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=wildcards.batch)

# prepare output adata
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_adata, compression='gzip')