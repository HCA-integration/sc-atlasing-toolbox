from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from scanpy.pp import pca
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    from rapids_singlecell.pp import harmony_integrate
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    from scanpy.external.pp import harmony_integrate
    logging.info('Importing rapids failed, using scanpy...')

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
keys = hyperparams.get('key', [])
if keys is None:
    keys = {batch_key}
elif isinstance(keys, str):
    keys = {batch_key, keys}
elif isinstance(keys, list):
    keys = set(keys).union({batch_key})
hyperparams['key'] = list(keys)

logging.info
(f'Read {input_file}...')
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
pca(adata, **pca_kwargs)
del adata.X
dask_compute(adata, layers=use_rep)

# run method
logging.info(f'Run harmonypy with parameters {pformat(hyperparams)}...')
harmony_integrate(
    adata,
    basis=use_rep,
    adjusted_basis='X_emb',
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