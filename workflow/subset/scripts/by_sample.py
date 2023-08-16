from utils import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
strategy = snakemake.wildcards.strategy
n_cell_max = snakemake.params.n_cells
sample_key = snakemake.params.sample_key
label_key = snakemake.params.label_key

print('read...')
adata = read_anndata(input_file)
print(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    print(adata)
    raise AssertionError(f'sample key "{sample_key}" not in adata')

# remove unannoted cells
adata = adata[adata.obs[label_key].notna()]
print('After removing unannotated cells:', adata.shape)

samples = []
n_cells = 0

shuffled_samples = adata.obs[sample_key].value_counts().sample(frac=1, random_state=42)
for sample, count in shuffled_samples.items():
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
print('write...')
adata.write(output_file)
