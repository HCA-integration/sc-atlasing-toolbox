from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from utils.processing import compute_neighbors
from utils.io import read_anndata, link_zarr


input_adata = snakemake.input[0]
output_file = snakemake.output[0]
# lineage_key = snakemake.wildcards.lineage_key
# lineage_specific = snakemake.wildcards.lineage_specific

files_to_overwrite = ['.zattrs', '.zgroup', 'obsp', 'var', 'uns']

logging.info('Read file...')
adata = read_anndata(input_adata)

# determine output types
# default output type is 'full'
output_type = adata.uns.get('integration', {'output_type': 'full'}).get('output_type')
output_types = [output_type] if isinstance(output_type, str) else output_type
adata.uns['output_types'] = output_types

for output_type in output_types:
    logging.info(f'Computing neighbors for output type {output_type}...')
    compute_neighbors(adata, output_type)
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
