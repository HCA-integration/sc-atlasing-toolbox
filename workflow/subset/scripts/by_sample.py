import anndata

input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
n_cell_max = snakemake.params.n_cells
strategy = snakemake.wildcards.strategy

adata = anndata.read_zarr(input_file)
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
adata.write_zarr(output_file)
