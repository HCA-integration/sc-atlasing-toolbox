from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, link_zarr
from utils.processing import compute_neighbors

input_file = snakemake.input[0]
output_file = snakemake.output[0]

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)
logging.info(adata.__str__())
logging.info(adata.uns)

output_type = adata.uns['integration']['output_type']
output_types = [output_type] if isinstance(output_type, str) else output_type

for output_type in output_types:
    logging.info(f'Computing neighbors for output type {output_type}...')
    compute_neighbors(adata, output_type)

logging.info('Write file...')
del adata.X
del adata.layers
adata.write_zarr(output_file)

input_files = [f.name for f in Path(input_file).iterdir()]
files_to_keep = [f for f in input_files if f not in ['obsm', 'obsp', 'uns', 'varm']]
link_zarr(
    in_dir=input_file,
    out_dir=output_file,
    file_names=files_to_keep,
    overwrite=True,
)
