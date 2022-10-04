import scib
import scanpy as sc

input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
method = snakemake.wildcards['method']
params = snakemake.params

adata = sc.read(input_adata)
adata_raw = adata.copy()

# run method
adata = scib.ig.scanorama(adata, batch=params['batch'])

# add metadata
adata.uns['dataset'] = params['dataset']

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

adata.layers['corrected_counts'] = adata.X.copy()
adata.layers['counts'] = adata_raw.layers['counts']
adata.layers['normcounts'] = adata_raw.layers['normcounts']

adata.write(output_adata, compression='gzip')