from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    from rapids_singlecell.pp import harmony_integrate
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    from scanpy.external.pp import harmony_integrate
    logging.info('Importing rapids failed, using scanpy...')

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
hyperparams = params.get('hyperparams')
if hyperparams is None:
    hyperparams = {}

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
logging.info(f'Run harmonypy with parameters {pformat(hyperparams)}...')
harmony_integrate(
    adata,
    key=wildcards.batch,
    use_rep=use_rep,
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