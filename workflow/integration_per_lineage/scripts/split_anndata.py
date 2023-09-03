import sys
from pathlib import Path
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.accessors import select_layer

input_file = snakemake.input[0]
output_dir = snakemake.output[0]

split_key = snakemake.params.lineage_key
batch_key = snakemake.params.batch
label_key = snakemake.params.label
# norm_layer = snakemake.params.norm_counts

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)

# remove unnecessary slots
del adata.obsm
del adata.obsp
del adata.varm
del adata.varp

# # ensure count matrix is correct
# logging.info('Select norm counts...')
# adata.X = select_layer(adata, norm_layer, dtype='float32')

# ensure count matrix is sparse
if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X, dtype='float32')


logging.info(f'Before filtering: {adata.shape}')
# remove splits with less than 100 cells
val_counts = adata.obs[split_key].value_counts()
adata = adata[adata.obs[split_key].isin(val_counts[val_counts > 100].index)]
adata = adata[adata.obs[split_key].notna()]

# remove unannotated cells
adata = adata[adata.obs[label_key].notna()]
logging.info(f'After filtering: {adata.shape}')

if adata.shape[0] == 0:
    raise AssertionError('No cells left after filtering')

splits = adata.obs[split_key].astype(str).unique()
logging.info(f'splits: {splits}')

highly_variable_genes_args = snakemake.params.get('hvg_args')
highly_variable_genes_args = {} if highly_variable_genes_args is None else highly_variable_genes_args
n_top_genes = highly_variable_genes_args.get('n_top_genes', 0)

for split in splits:

    logging.info(f'Split by {split_key}={split}')
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split]

    # # ensure enough cells per batch
    # n_cells_before = adata_sub.n_obs
    # logging.info(f'number of cells before filtering: {n_cells_before}')

    # val_counts = adata_sub.obs[batch_key].value_counts()
    # batches_to_keep = val_counts[val_counts > n_top_genes].index
    # adata_sub = adata_sub[adata_sub.obs[batch_key].isin(batches_to_keep)]
    # logging.info(f'number of cells after filtering: {adata_sub.n_obs}, removed {adata_sub.n_obs-n_cells_before} cells')

    if adata_sub.n_obs == 0:
        logging.info('No cells left after filtering batches by HVG, skipping...')
        continue

    # write to file
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"lineage~{split_file}.zarr"

    logging.info(f'write to {out_file}...')
    adata_sub.write_zarr(out_file)
    del adata_sub
