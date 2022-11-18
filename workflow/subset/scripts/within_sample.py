import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
n_cell_max = snakemake.params.n_cells
strategy = snakemake.wildcards.strategy

adata = sc.read(input_h5ad)
print(f'Shape before filtering: {adata.shape}')

adatas = []
n_samples = adata.obs['sample'].nunique()
n_cells_per_sample = int(n_cell_max / n_samples)

print(f'Subset to {n_cells_per_sample} cells per sample...')

for sample in adata.obs['sample'].unique():
    ad = adata[adata.obs['sample'] == sample]
    sc.pp.subsample(ad, n_obs=n_cells_per_sample)
    adatas.append(ad)

adata = sc.concat(adatas)
adata.uns['subset'] = strategy
print(f'Subset to {adata.n_obs} cells, {n_samples} samples')
print(adata)

# save
adata.write(output_h5ad)