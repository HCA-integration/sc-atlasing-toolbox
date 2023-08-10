"""
Highly variable gene selection
- lineage specific HVGs
"""
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params['args']
batch_key = snakemake.params['batch']
lineage_key = snakemake.params['lineage']

if args is None:
    args = {}
logging.info(str(args))

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

adata.uns["log1p"] = {"base": None}
sc.pp.filter_genes(adata, min_cells=1)

if args is False:
    logging.info('No highly variable gene parameters provided, including all genes...')
    adata.var['highly_variable'] = True
else:
    if lineage_key is not None:
        logging.info(f'lineage-specific highly variable gene selection using "{lineage_key}"')
        if batch_key in adata.obs.columns:
            adata.obs['hvg_batch'] = adata.obs[batch_key].astype(str) + '_' + adata.obs[lineage_key].astype(str)
        else:
            adata.obs['hvg_batch'] = adata.obs[lineage_key]
        batch_key = 'hvg_batch'

    logging.info('Select features...')
    sc.pp.highly_variable_genes(adata, batch_key=batch_key, **args)

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['highly_variable_genes'] = args

# remove counts
del adata.X
del adata.layers

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)
