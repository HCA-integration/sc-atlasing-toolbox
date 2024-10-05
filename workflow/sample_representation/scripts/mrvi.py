import logging
import warnings
import torch

import patient_representation as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
use_rep = snakemake.params.get('use_rep')
n_epochs = snakemake.params.get('n_epochs')
n_epochs = int(n_epochs) if n_epochs is not None else n_epochs

use_gpu = torch.cuda.is_available()
print(f'GPU available: {use_gpu}', flush=True)

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_zarr,
    X=use_rep,
    obs='obs',
    var='var',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

logging.info(f'Calculating MrVI representation for "{cell_type_key}", using cell features from "{use_rep}"')
representation_method = pr.tl.MrVI(
    sample_key=sample_key,
    cells_type_key=cell_type_key,
    layer='X',
    max_epochs=n_epochs,
    accelerator='gpu' if use_gpu else 'auto',
)
representation_method.prepare_anndata(adata)
distances = representation_method.calculate_distance_matrix(force=True)

# Create empty dataframe to put sample names to obs_names and var_names
samples_df = pd.DataFrame(index=representation_method.samples)

output_adata = sc.AnnData(X=distances, var=samples_df, obs=samples_df)

logging.info(output_adata.__str__())

logging.info(f'Write "{output_zarr}"...')
output_adata.write_zarr(output_zarr)
