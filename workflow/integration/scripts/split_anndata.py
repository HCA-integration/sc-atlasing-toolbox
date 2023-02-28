import sys
from pathlib import Path
from scipy.sparse import csr_matrix
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

from methods.utils import read_anndata

input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.lineage_key
batch_key = snakemake.wildcards.batch

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

adata = read_anndata(input_file)

# preprocessing
if 'preprocessing' in adata.uns:
    n_top_genes = adata.uns['preprocessing']['hvg']
else:
    n_top_genes = 2000
    adata.uns['preprocessing'] = {'hvg': n_top_genes}

splits = adata.obs[split_key].unique()
print(f'splits: {splits}', file=sys.stderr)

for split in splits:
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split].copy()
    adata_sub.uns["log1p"] = {"base": None}
    adata_sub.X = csr_matrix(adata_sub.X)

    if n_top_genes > 0:
        sc.pp.filter_genes(adata_sub, min_cells=1)
        print(f'HVGs to {n_top_genes} genes...', file=sys.stderr)
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_top_genes, batch_key=batch_key)
        print('PCA...', file=sys.stderr)
        sc.pp.pca(adata_sub, use_highly_variable=True, svd_solver='arpack')
    else:
        sc.pp.pca(adata_sub, use_highly_variable=False, svd_solver='arpack')

    print('compute neighbors...', file=sys.stderr)
    sc.pp.neighbors(adata_sub, use_rep='X_pca')

    # write to file
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"{split_file}.h5ad"

    print(f'write to {out_file}...', file=sys.stderr)
    adata_sub.write(out_file, compression='lzf')
