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
import anndata as ad
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
import sparse
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError:
    sc = scanpy
    logging.info('Importing rapids failed, using scanpy...')
    rapids = False

from utils.io import read_anndata, write_zarr_linked, csr_matrix_int64_indptr
from utils.misc import apply_layers


def filter_genes(adata):
    logging.info('Filter genes...')
    if isinstance(adata.X, da.Array):
        adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
    gene_subset, _ = scanpy.pp.filter_genes(adata.X, min_cells=1)
    if isinstance(gene_subset, da.Array):
        gene_subset = gene_subset.compute()
    if any(gene_subset == False):
        logging.info(f'Subset to {sum(gene_subset)}/{adata.n_vars} filtered genes...')
        adata = adata[:, gene_subset].copy()
    if isinstance(adata.X, da.Array):
        adata.X = adata.X.map_blocks(csr_matrix_int64_indptr).compute()
    return adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
batch_key = snakemake.params.get('batch')
lineage_key = snakemake.params.get('lineage')
dask = snakemake.params.get('dask', False) and not rapids
backed = snakemake.params.get('backed', False) and dask and not rapids

if args is None:
    args = {}
subset_to_hvg = isinstance(args, dict) and args.get('subset', False)
logging.info(str(args))

logging.info(f'Read {input_file}...')
kwargs = dict(X='X', obs='obs', var='var', uns='uns', backed=backed, dask=dask)
if subset_to_hvg:
    kwargs |= dict(layers='layers')
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())
var = adata.var.copy()

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
    var['highly_variable'] = True
else:
    if lineage_key is not None:
        logging.info(f'lineage-specific highly variable gene selection using "{lineage_key}"')
        if batch_key in adata.obs.columns: # combining lineage and batch
            adata.obs['hvg_batch'] = adata.obs[batch_key].astype(str) + '_' + adata.obs[lineage_key].astype(str)
        else: # only taking lineage if no batch key is specified
            adata.obs['hvg_batch'] = adata.obs[lineage_key]
        batch_key = 'hvg_batch'

    adata = filter_genes(adata)
    
    # make sure data is on GPU for rapids_singlecell
    if rapids:
        sc.utils.anndata_to_GPU(adata)
    
    logging.info(f'Select features with arguments: {args}...')
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
    for column, default_value in hvg_column_map.items():
        if column not in adata.var.columns:
            continue
        dtype = adata.var[column].dtype
        var[column] = default_value
        var[column] = var[column].astype(dtype)
        var.loc[adata.var_names, column] = adata.var[column]

logging.info(f'Write to {output_file}...')
files_to_keep = ['uns', 'var']
if subset_to_hvg:
    logging.info('Subset to highly variable genes...')
    adata = adata[:, adata.var['highly_variable']].copy()
    files_to_keep.extend(['X', 'layers', 'varm', 'varp'])
    logging.info(adata.__str__())
else:
    adata = ad.AnnData(var=var, uns=adata.uns)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep
)