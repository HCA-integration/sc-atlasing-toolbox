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

samples = []
n_cells = 0

for sample, count in adata.obs[sample_key].value_counts(sort=False).items():
    n_cells += count
    print('sample:', sample)
    print('count:', count)
    if n_cells > n_cell_max:
        break
    samples.append(sample)

print(f'Subset to {n_cells} cells, {len(samples)} samples...')

adata = adata[adata.obs[sample_key].isin(samples)]
adata.uns['subset'] = strategy
print(adata)

# save
adata.write_zarr(output_file)
