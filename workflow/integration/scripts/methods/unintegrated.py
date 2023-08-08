import scanpy as sc

from utils import add_metadata, read_anndata, process, select_layer

input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = read_anndata(input_adata)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# prepare output adata
adata = adata_raw
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
adata.obsm['X_emb'] = adata.obsm['X_pca']
sc.pp.neighbors(adata)
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)
