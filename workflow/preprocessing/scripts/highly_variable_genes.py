"""
Highly variable gene selection
- lineage specific HVGs
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
from scipy import sparse
import scanpy as sc
from utils.io import read_anndata, link_zarr


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
batch_key = snakemake.params.get('batch')
lineage_key = snakemake.params.get('lineage')

if args is None:
    args = {}
logging.info(str(args))

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

adata.uns["log1p"] = {"base": None}

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
    adata_hvg = adata.copy()
    sc.pp.filter_genes(adata_hvg, min_cells=1)
    sc.pp.highly_variable_genes(adata_hvg, batch_key=batch_key, **args)

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
        if column in adata_hvg.var.columns
    ]
    adata.var[hvg_columns] = adata_hvg.var[hvg_columns]
    adata.var = adata.var.fillna(hvg_column_map)

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['highly_variable_genes'] = args

# remove counts
del adata.X
del adata.layers

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    files_to_keep = ['uns', 'var']
    if isinstance(args, dict) and 'subset' in args and args['subset']:
        logging.info('Data subsetted to highly variable genes, keep matrices, varm and varp...')
        files_to_keep.extend(['layers', 'X', 'varm', 'varp'])

    input_files = [f.name for f in Path(input_file).iterdir()]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=[f for f in input_files if f not in files_to_keep],
        overwrite=True,
)
