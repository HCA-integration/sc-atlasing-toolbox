import scanpy as sc

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
n_cell_max = snakemake.params.n_cells
strategy = snakemake.wildcards.strategy

adata = sc.read(input_h5ad)
print(f'Shape before filtering: {adata.shape}')

samples = []
n_cells = 0

for sample, count in adata.obs['sample'].value_counts(sort=False).items():
    n_cells += count
    samples.append(sample)
    if n_cells > n_cell_max:
        break

print(f'Subset to {n_cells} cells, {len(samples)} samples...')

adata = adata[adata.obs['sample'].isin(samples)]
adata.uns['subset'] = strategy
print(adata)

# save
adata.write(output_h5ad)