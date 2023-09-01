import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import scanpy as sc

from utils import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
strategy = snakemake.wildcards.strategy
n_cell_max = snakemake.params.n_cells
sample_key = snakemake.params.sample_key


logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)
logging.info(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    logging.info(adata.__str__())
    raise AssertionError(f'sample key "{sample_key}" not in adata')

adatas = []
n_samples = adata.obs[sample_key].nunique()
n_cells_per_sample = int(n_cell_max / n_samples)

logging.info(f'Subset to {n_cells_per_sample} cells per sample...')

for sample in adata.obs[sample_key].unique():
    ad = adata[adata.obs[sample_key] == sample]
    sc.pp.subsample(ad, n_obs=np.min([n_cells_per_sample, ad.n_obs]))
    adatas.append(ad)

adata = sc.concat(adatas)
adata.uns['subset'] = strategy
logging.info(f'Subset to {adata.n_obs} cells, {n_samples} samples')
logging.info(adata.__str__())

# save
logging.info('Write...')
adata.write_zarr(output_file)
