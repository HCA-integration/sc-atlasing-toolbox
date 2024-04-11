from pathlib import Path
from scipy.sparse import csr_matrix
import numpy as np
import scanpy as sc
from anndata import AnnData

adata = sc.datasets.pbmc68k_reduced()
adata.X = adata.raw.X.copy().todense()

np.random.seed(42)
n_batch = 3
adata.obs['batch'] = np.random.randint(0, n_batch, adata.n_obs)
# add batch effect to counts
for i in range(n_batch):
    adata[adata.obs['batch'] == i].X += i
adata.obs['batch'] = adata.obs['batch'].astype(str)
adata.obs['batch_2'] = adata.obs['phase']

adata.layers['counts'] = csr_matrix(np.exp(adata.X).astype(int) - 1)
adata.X = csr_matrix(adata.X)
adata.layers['normcounts'] = adata.X.copy()

# add column with NA columns
adata.obs['na_column'] = adata.obs['bulk_labels'].astype(str)
adata.obs.loc[adata.obs['na_column'] == 'Dendritic', 'na_column'] = np.nan
adata.obs['na_str_column'] = adata.obs['na_column'].astype(str)

# add bolean columns
adata.obs['is_cd14_mono'] = adata.obs['bulk_labels'] == 'CD14+ Monocyte'

# add var mask
adata.var['highly_variable_2'] = adata.var['highly_variable'].copy()

# add raw
adata.raw = AnnData(
    X=adata.layers['counts'],
    obs=adata.obs,
    var=adata.var,
)

print(adata)

out_file = Path(__file__).parent / 'pbmc68k.h5ad'
print(f'writing to {out_file}...')
adata.write(out_file, compression='gzip')
