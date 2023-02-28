from pathlib import Path
from scipy.sparse import csr_matrix
import numpy as np
import scanpy as sc

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

adata.layers['counts'] = csr_matrix(np.exp(adata.X) - 1)
adata.X = csr_matrix(adata.X)
adata.layers['normcounts'] = adata.X.copy()
del adata.raw
del adata.uns

print(adata)

out_file = Path(__file__).parent / 'pbmc68k.h5ad'
print(f'writing to {out_file}...')
adata.write(out_file, compression='gzip')
