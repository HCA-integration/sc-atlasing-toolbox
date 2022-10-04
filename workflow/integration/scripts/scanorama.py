import scib
import scanpy as sc

input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
method = snakemake.wildcards['method']
params = snakemake.params

adata = sc.read(input_adata)

# add metadata
if 'methods' in adata.uns:
    adata.uns['methods'].append(method)
else:
    adata.uns['methods'] = [method]

adata.uns['integration'] = {
    'method': method,
    'label_key': params['label'],
    'batch_key': params['batch'],
    'output_type': params['output_type']
}

# run method
scib.ig.scanorama(adata, batch=params['batch'])

adata.write(output_adata, compression='gzip')