from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
hyperparams = params.get('hyperparams')
if hyperparams is None:
    hyperparams = {}
if 'covariates' in hyperparams:
    covariates = hyperparams['covariates']
    if isinstance(covariates, str):
        hyperparams['covariates'] = [covariates]

adata = read_anndata(
    input_file,
    X='layers/norm_counts',
    obs='obs',
    var='var',
    uns='uns'
)

# run method
logging.info(f'Run Combat with parameters {pformat(hyperparams)}...')
sc.pp.combat(
    adata,
    key=wildcards.batch,
    **hyperparams
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['X', 'uns'])