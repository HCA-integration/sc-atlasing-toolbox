from pathlib import Path
import logging
logger = logging.getLogger('Preprocess for metrics')
logger.setLevel(logging.INFO)

from utils.processing import compute_neighbors
from utils.io import read_anndata, link_zarr


input_adata = snakemake.input[0]
output_file = snakemake.output[0]
label_key = snakemake.params.label_key
neighbor_args = snakemake.params.neighbor_args

files_to_overwrite = ['obsp', 'var', 'uns']

logging.info('Read file...')
adata = read_anndata(input_adata)  # TODO: read only necessary slots

# determine output types
# default output type is 'full'
output_type = adata.uns.get('output_type', 'full')

n_obs = adata.n_obs
# logger.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# TODO: Error when connectivities/distances have different dimensions than X
# logger.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logger.info(f'After: {adata.n_obs} cells')
force_neighbors = n_obs > adata.n_obs

logging.info(f'Computing neighbors for output type {output_type}...')
compute_neighbors(adata, output_type, force=force_neighbors, **neighbor_args)
if output_type == 'full':
    files_to_overwrite.extend(['obsm', 'varm'])

# write to file
logging.info(adata.__str__())
adata.write_zarr(output_file)

if input_adata.endswith('.zarr'):
    input_file = Path(input_adata)
    input_zarr_files = [f.name for f in Path(input_file).iterdir()]
    files_to_link = [f for f in input_zarr_files if f not in files_to_overwrite]
    link_zarr(
        in_dir=input_file,
        out_dir=Path(output_file),
        file_names=files_to_link,
        overwrite=True,
)
