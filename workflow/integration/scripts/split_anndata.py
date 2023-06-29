import sys
from pathlib import Path
from scipy.sparse import csr_matrix
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from methods.utils import read_anndata

input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.lineage_key
batch_key = snakemake.wildcards.batch
label_key = snakemake.params.label

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)

# remove unannoted cells
logging.info(f'Before filtering: {adata.shape}')
val_counts = adata.obs[split_key].value_counts()
adata = adata[adata.obs[split_key].isin(val_counts[val_counts > 100].index)]
adata = adata[adata.obs[split_key].notna()]
adata = adata[adata.obs[label_key].notna()]
logging.info(f'After filtering: {adata.shape}')

if adata.shape[0] == 0:
    raise AssertionError('No cells left after filtering')

splits = adata.obs[split_key].astype(str).unique()
logging.info(f'splits: {splits}')

# preprocessing
if 'preprocessing' in adata.uns:
    highly_variable_genes_args = adata.uns['preprocessing']['highly_variable_genes']

n_top_genes = 0
if 'highly_variable_genes' in adata.var.columns:
    n_top_genes = adata.var['highly_variable'].sum()

for split in splits:
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split].copy()

    # preprocessing
    adata_sub.uns["log1p"] = {"base": None}
    adata_sub.X = csr_matrix(adata_sub.X)

    if n_top_genes > 0:
        sc.pp.filter_genes(adata_sub, min_cells=1)
        logging.info(f'HVGs to {n_top_genes} genes...')
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_top_genes, batch_key=batch_key)
        logging.info('PCA...')
        sc.pp.pca(adata_sub, use_highly_variable=True, svd_solver='arpack')
    else:
        sc.pp.pca(adata_sub, use_highly_variable=False, svd_solver='arpack')

    logging.info('compute neighbors...')
    sc.pp.neighbors(adata_sub, use_rep='X_pca')

    # write to file
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"{split_file}.h5ad"

    logging.info(f'write to {out_file}...')
    adata_sub.write(out_file, compression='lzf')
