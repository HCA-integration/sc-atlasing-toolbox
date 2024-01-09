import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

files_to_keep = ['obsm', 'obsp', 'uns']

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    obs='obs',
    var='var',
    obsm='obsm',
    uns='uns'
)

assert 'X_pca' in adata.obsm.keys(), 'PCA is missing'

# quickfix: remove batches with fewer than 3 cells
min_batches = adata.obs.groupby(batch_key).filter(lambda x: len(x) > 3).index
if min_batches.nunique() < adata.n_obs:
    files_to_keep.extend(['obs'])
    # adata.layers = read_anndata(input_file, layers='layers').layers
    adata = adata[min_batches]

# run method
adata = sc.external.pp.bbknn(adata, batch_key=batch_key, use_rep='X_pca', copy=True)

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