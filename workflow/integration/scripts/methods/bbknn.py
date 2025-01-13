from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc

from integration_utils import clean_categorical_column, add_metadata, \
    remove_slots, get_hyperparams, PCA_PARAMS
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams

pca_kwargs, hyperparams = get_hyperparams(
    hyperparams=hyperparams,
    model_params=PCA_PARAMS,
)
hyperparams = {'pynndescent_random_state': params.get('seed', 0)} | hyperparams

files_to_keep = ['obsm', 'obsp', 'uns']

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/normcounts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

# subset features
adata, _ = subset_hvg(
    adata,
    var_column='integration_features',
    compute_dask=True
)

# recompute PCA according to user-defined hyperparameters
logging.info(f'Compute PCA with parameters {pformat(pca_kwargs)}...')
use_rep = 'X_pca'
# adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
sc.pp.pca(adata, **pca_kwargs)
del adata.X
# dask_compute(adata, layers=use_rep)

# quickfix: remove batches with fewer than 3 cells
neighbors_within_batch = hyperparams.get('neighbors_within_batch', 3)
min_batches = adata.obs.groupby(
    batch_key,
    observed=False
).filter(
    lambda x: len(x) > neighbors_within_batch
).index

if min_batches.nunique() < adata.n_obs:
    files_to_keep.extend(['obs'])
    # adata.layers = read_anndata(input_file, layers='layers').layers
    adata = adata[min_batches]

# run method
logging.info(f'Run BBKNN with parameters {pformat(hyperparams)}...')
sc.external.pp.bbknn(
    adata,
    batch_key=batch_key,
    use_rep=use_rep,
    copy=False,
    **hyperparams,
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
    files_to_keep=files_to_keep,
)