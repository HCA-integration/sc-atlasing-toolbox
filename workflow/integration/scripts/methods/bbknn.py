import scanpy as sc

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

adata = read_anndata(input_adata)
adata.X = select_layer(adata, params['norm_counts'])

# subset to HVGs
adata = adata[:, adata.var['highly_variable']]

# quickfix: remove batches with fewer than 3 cells
adata = adata[adata.obs.groupby(batch_key).filter(lambda x: len(x) > 3).index]

# run method
adata = sc.external.pp.bbknn(adata, batch_key=batch_key, use_rep='X_pca', copy=True)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)
