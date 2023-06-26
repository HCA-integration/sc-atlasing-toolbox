import scanpy as sc

from utils import add_metadata, read_anndata, process, select_layer


input_adata = snakemake.input.h5ad
output_adata = snakemake.output.h5ad
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

adata_raw = read_anndata(input_adata)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# quickfix: remove batches with fewer than 3 cells
adata_raw = adata_raw[adata_raw.obs.groupby(batch_key).filter(lambda x: len(x) > 3).index]

# run method
adata = sc.external.pp.bbknn(adata_raw, batch_key=batch_key, use_rep='X_pca', copy=True)

# prepare output adata
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_adata, compression='gzip')
