import sys
from pathlib import Path
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from methods.utils import read_anndata, select_layer

input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]

split_key = snakemake.wildcards.lineage_key
batch_key = snakemake.wildcards.batch
label_key = snakemake.params.label
norm_layer = snakemake.params.norm_counts

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)

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

# get preprocessing parameters
try:
    highly_variable_genes_args = adata.uns['preprocessing']['highly_variable_genes']
except KeyError:
    highly_variable_genes_args = {}
n_top_genes = adata.var['highly_variable'].sum() if 'highly_variable' in adata.var.columns else None
highly_variable_genes_args.update(
    dict(n_top_genes=n_top_genes, batch_key=batch_key)
)

for split in splits:

    logging.info(f'Split by {split_key}={split}')
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split].copy()

    # ensure enough cells per batch
    val_counts = adata_sub.obs[batch_key].value_counts()
    adata_sub = adata_sub[adata_sub.obs[batch_key].isin(val_counts[val_counts > n_top_genes].index)]

    # ensure count matrix is correct
    #logging.info('Select norm counts...')
    #adata_sub.X = select_layer(adata_sub, norm_layer, force_dense=True, dtype='float32')

    # preprocessing
    adata_sub.uns["log1p"] = {"base": None}

    if n_top_genes:
        sc.pp.filter_genes(adata_sub, min_cells=1)

        logging.info(f'HVGs to {n_top_genes} genes...')
        sc.pp.highly_variable_genes(adata_sub, **highly_variable_genes_args)

        logging.info('PCA...')
        sc.pp.pca(
            adata_sub,
            use_highly_variable=True,
            svd_solver='arpack'
        )
    else:
        logging.info('PCA...')
        sc.pp.pca(
            adata_sub,
            use_highly_variable=False,
            svd_solver='arpack'
        )

    logging.info('Compute neighbors...')
    try:
        sc.pp.neighbors(adata_sub, method='rapids', use_rep='X_pca')
    except Exception as e:
        logging.info(e)
        logging.info('Rapids failed, defaulting to UMAP implementation')
        sc.pp.neighbors(adata_sub, use_rep='X_pca')

    # write to file
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"{split_file}.h5ad"

    logging.info(f'write to {out_file}...')
    adata_sub.X = sparse.csr_matrix(adata_sub.X)
    adata_sub.write(out_file, compression='lzf')
    del adata_sub
