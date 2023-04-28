"""
Highly variable gene selection
- lineage specific HVGs
"""
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
print(args)

print('read...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    adata.write(output_file)
    exit(0)

adata.uns["log1p"] = {"base": None}
sc.pp.filter_genes(adata, min_cells=1)

if lineage_key is None:
    sc.pp.highly_variable_genes(adata, batch_key=batch_key, **args)
else:
    print(f'lineage-specific highly variable gene selection using "{lineage_key}"')
    if batch_key in adata.obs.columns:
        adata.obs['hvg_batch'] = adata.obs[batch_key].astype(str) + '_' + adata.obs[lineage_key].astype(str)
    else:
        adata.obs['hvg_batch'] = adata.obs[lineage_key]
    sc.pp.highly_variable_genes(adata, batch_key='hvg_batch', **args)

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['highly_variable_genes'] = args

print('write...')
adata.X = sparse.csr_matrix(adata.X)
adata.write(output_file, compression='lzf')
