import logging
import warnings

import patient_representation as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

input_file = snakemake.input.zarr
prepare_file = snakemake.input.prepare
output_file = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
use_rep = snakemake.params.get('use_rep')
var_mask = snakemake.params.get('var_mask')

logging.info(f'Read "{input_file}"...')
n_obs = read_anndata(input_file, obs='obs').n_obs
dask = n_obs > 1e6
adata = read_anndata(
    input_file,
    X=use_rep,
    obs='obs',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

# subset HVGs
if var_mask is not None:
    adata.var = read_anndata(input_file, var='var').var
    adata = adata[:, adata.var[var_mask]].copy()
dask_compute(adata)

logging.info(f'Calculating TotalPseudobulk representation for sample_key={sample_key}, using cell features from "{use_rep}"')
representation_method = pr.tl.TotalPseudobulk(
    sample_key=sample_key,
    cells_type_key=sample_key,
    layer='X',
)
# representation_method.prepare_anndata(adata)
representation_method.adata = adata
representation_method.samples = representation_method.adata.obs[sample_key].unique()

# create new AnnData object for patient representations
adata = sc.AnnData(obs=pd.DataFrame(index=representation_method.samples))
adata.obsm['distances'] = representation_method.calculate_distance_matrix(force=True)
adata.obsm['X_emb'] = representation_method.patient_representations
samples = read_anndata(prepare_file, obs='obs').obs_names
adata = adata[samples].copy()

# compute kNN graph
# sc.pp.neighbors(adata, use_rep='distances', metric='precomputed', transformer='sklearn')
# sc.pp.neighbors(adata, use_rep='X_emb', key_added='X_emb')
sc.pp.neighbors(adata, use_rep='X_emb')

logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=prepare_file,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns']
)