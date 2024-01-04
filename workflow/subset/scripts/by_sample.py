import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
strategy = snakemake.wildcards.strategy
n_cell_max = snakemake.params.get('n_cells')
n_cell_max = np.iinfo(int).max if n_cell_max is None else int(n_cell_max)
sample_key = snakemake.params.get('sample_key')
# label_key = snakemake.params.get('label_key')


logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, backed=True)
logging.info(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    logging.info(adata.__str__())
    raise AssertionError(f'sample key "{sample_key}" not in adata')

# remove unannoted cells
# adata = adata[adata.obs[label_key].notna()]
# logging.info(f'After removing unannotated cells: {adata.shape}')

samples = []
n_cells = 0

shuffled_samples = adata.obs[sample_key].value_counts().sample(frac=1, random_state=42)
for sample, count in shuffled_samples.items():
    n_cells += count
    logging.info(f'sample: {sample}')
    logging.info(f'count: {count}')
    if len(samples) > 0 and n_cells > n_cell_max:
        break
    samples.append(sample)

logging.info(f'Subset to {n_cells} cells, {len(samples)} samples...')

adata = adata[adata.obs[sample_key].isin(samples)]
adata.uns['subset'] = strategy
logging.info(adata)

# save
logging.info('Write...')
adata.write_zarr(output_file)
