from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import torch
from harmony import harmonize

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
hyperparams = params.get('hyperparams')
if hyperparams is None:
    hyperparams = {}

# check GPU
use_gpu = torch.cuda.is_available()
logging.info(f'GPU available: {use_gpu}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    obs='obs',
    var='var',
    obsm='obsm',
    uns='uns'
)

use_rep = hyperparams.pop('use_rep', 'X_pca')
assert use_rep in adata.obsm.keys(), f'{use_rep} is missing'

# run method
logging.info(f'Run Harmony pytorch with parameters {pformat(hyperparams)}...')
adata.obsm['X_emb'] = harmonize(
    X=adata.obsm[use_rep],
    batch_mat=adata.obs,
    batch_key=wildcards.batch,
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