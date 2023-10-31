import scib

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata = read_anndata(input_adata)
adata.X = select_layer(adata, params['norm_counts'])

# subset to HVGs
adata = adata[:, adata.var['highly_variable']]

# run method
adata = scib.ig.scanorama(adata, batch=wildcards.batch)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)