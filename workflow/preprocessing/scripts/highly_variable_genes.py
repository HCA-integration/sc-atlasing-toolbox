"""
Highly variable gene selection
- lineage specific HVGs
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
import scanpy
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError:
    sc = scanpy
    logging.info('Importing rapids failed, using scanpy...')
    rapids = False

from utils.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
batch_key = snakemake.params.get('batch')
lineage_key = snakemake.params.get('lineage')

if args is None:
    args = {}
subset_to_hvg = isinstance(args, dict) and args.get('subset', False)
logging.info(str(args))

logging.info(f'Read {input_file}...')
kwargs = dict(X='X', obs='obs', var='var', uns='uns')
if subset_to_hvg:
    kwargs |= dict(layers='layers')
adata = read_anndata(input_file, **kwargs)
var = adata.var

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['highly_variable_genes'] = args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.var['highly_variable'] = True
    adata.write_zarr(output_file)
    exit(0)

# scanpy.pp.log1p was supposed to add it but it's not saved
adata.uns["log1p"] = {"base": None}

if args is False:
    logging.info('No highly variable gene parameters provided, including all genes...')
    adata.var['highly_variable'] = True
else:
    if lineage_key is not None:
        logging.info(f'lineage-specific highly variable gene selection using "{lineage_key}"')
        if batch_key in adata.obs.columns: # combining lineage and batch
            adata.obs['hvg_batch'] = adata.obs[batch_key].astype(str) + '_' + adata.obs[lineage_key].astype(str)
        else: # only taking lineage if no batch key is specified
            adata.obs['hvg_batch'] = adata.obs[lineage_key]
        batch_key = 'hvg_batch'

    logging.info('Filter genes...')
    gene_subset, _ = scanpy.pp.filter_genes(adata, min_cells=1, inplace=False)
    if any(gene_subset == False):
        logging.info('Subset to filtered genes...')
        adata = adata[:, gene_subset]
    
    # make sure data is on GPU for rapids_singlecell
    if rapids:
        sc.utils.anndata_to_GPU(adata)
    
    logging.info('Select features...')
    if 'subset' in args:
        del args['subset']
    sc.pp.highly_variable_genes(
        adata,
        batch_key=batch_key,
        **args
    )

    # add HVG info back to adata
    hvg_column_map = {
        'highly_variable': False,
        'means': 0,
        'dispersions': 0,
        'dispersions_norm': 0,
        'highly_variable_nbatches': 0,
        'highly_variable_intersection': False,
    }
    hvg_columns = [
        column
        for column in hvg_column_map
        if column in adata.var.columns
    ]
    var[hvg_columns] = adata.var[hvg_columns]
    var = var.fillna(hvg_column_map)

logging.info(f'Write to {output_file}...')
files_to_keep = ['uns', 'var']
if subset_to_hvg:
    logging.info('Subset to highly variable genes...')
    adata = adata[:, adata.var['highly_variable']].copy()
    files_to_keep.extend(['X', 'layers', 'varm', 'varp'])
else:
    del adata.X
    adata.var = var

logging.info(adata.__str__())

write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep
)