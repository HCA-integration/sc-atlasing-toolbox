import scib
import scanpy as sc

from utils import process


input_adata = snakemake.input.h5ad
output_adata = snakemake.output.h5ad
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = sc.read(input_adata)
adata_raw.X = adata_raw.layers['normcounts'].copy()

# run method
adata = scib.ig.bbknn(adata_raw, batch=wildcards.batch)
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])

# add metadata
adata.uns['dataset'] = wildcards.dataset

if 'methods' in adata.uns:
    adata.uns['methods'].append(wildcards.method)
else:
    adata.uns['methods'] = [wildcards.method]

adata.uns['integration'] = {
    'method': wildcards.method,
    'label_key': wildcards.label,
    'batch_key': wildcards.batch,
    'output_type': params['output_type']
}

adata.write(output_adata, compression='gzip')
