import numpy as np
import scanpy as sc

from utils import read_anndata

input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
strategy = snakemake.wildcards.strategy
n_cell_max = snakemake.params.n_cells
sample_key = snakemake.params.sample_key

adata = read_anndata(input_file)
print(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    print(adata)
    raise AssertionError(f'sample key "{sample_key}" not in adata')

adatas = []
n_samples = adata.obs[sample_key].nunique()
n_cells_per_sample = int(n_cell_max / n_samples)

print(f'Subset to {n_cells_per_sample} cells per sample...')

for sample in adata.obs[sample_key].unique():
    ad = adata[adata.obs[sample_key] == sample]
    sc.pp.subsample(ad, n_obs=np.min([n_cells_per_sample, ad.n_obs]))
    adatas.append(ad)

adata = sc.concat(adatas)
adata.uns['subset'] = strategy
print(f'Subset to {adata.n_obs} cells, {n_samples} samples')
print(adata)

# save
adata.write_zarr(output_file)
