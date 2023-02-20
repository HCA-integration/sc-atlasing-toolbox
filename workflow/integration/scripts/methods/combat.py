import scib
import scanpy as sc

from utils import add_metadata, read_anndata, process


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = read_anndata(input_adata)
adata_raw.X = adata_raw.layers['normcounts'].copy()

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# run method
adata = scib.ig.combat(adata_raw, batch=wildcards.batch)

# prepare output adata
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_adata, compression='gzip')