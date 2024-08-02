from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import torch
import scanpy as sc
from harmony import harmonize

from integration_utils import add_metadata, remove_slots, get_hyperparams, PCA_PARAMS
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
var_mask = wildcards.var_mask

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams

pca_kwargs, hyperparams = get_hyperparams(
    hyperparams=hyperparams,
    model_params=PCA_PARAMS,
)
hyperparams = {'random_state': params.get('seed', 0)} | hyperparams


# set harmony var_use
keys = hyperparams.get('batch_key', [])
if keys is None:
    keys = {batch_key}
elif isinstance(keys, str):
    keys = {batch_key, keys}
elif isinstance(keys, list):
    keys = set(keys).union({batch_key})
hyperparams['batch_key'] = list(keys)

# check GPU
use_gpu = torch.cuda.is_available()
logging.info(f'GPU available: {use_gpu}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/norm_counts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

# subset features
adata, subsetted = subset_hvg(adata, var_column=var_mask, compute_dask=False)

# recompute PCA according to user-defined hyperparameters
logging.info(f'Compute PCA with parameters {pformat(pca_kwargs)}...')
use_rep = 'X_pca'
adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
sc.pp.pca(adata, **pca_kwargs)
del adata.X
dask_compute(adata, layers=use_rep)

# run method
logging.info(f'Run Harmony pytorch with parameters {pformat(hyperparams)}...')
adata.obsm['X_emb'] = harmonize(
    X=adata.obsm[use_rep],
    batch_mat=adata.obs,
    use_gpu=use_gpu,
    n_jobs=snakemake.threads,
    **hyperparams
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)