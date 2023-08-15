from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from metrics.utils import anndata_to_mudata, compute_neighbors, get_from_adata
from utils.io import read_anndata_or_mudata, link_zarr


input_adata = snakemake.input[0]
output_file = snakemake.output[0]
lineage_key = snakemake.wildcards.lineage_key
lineage_specific = snakemake.wildcards.lineage_specific

files_to_overwrite = ['.zattrs', '.zgroup', 'obsp', 'var', 'uns']

logging.info('Read file...')
adata = read_anndata_or_mudata(input_adata)
meta = get_from_adata(adata)

mudata = anndata_to_mudata(adata, group_key=lineage_key)

for lineage in mudata.mod:
    logging.info(f'Processing lineage {lineage}...')
    ad = mudata[lineage]
    for output_type in meta['output_types']:
        logging.info(f'Computing neighbors for output type {output_type}...')
        compute_neighbors(ad, output_type)
        if output_type == 'full':
            files_to_overwrite.extend(['obsm', 'varm'])

# write to file
logging.info(mudata.__str__())
mudata.write_zarr(output_file)

if lineage_specific == 'global':
    exit()

for lineage in mudata.mod:
    input_file = Path(input_adata)
    if input_adata.endswith('mu.zarr'):
        input_file = input_file / 'mod' / lineage
    input_zarr_files = [f.name for f in Path(input_file).iterdir()]
    files_to_link = [f for f in input_zarr_files if f not in files_to_overwrite]
    link_zarr(
        in_dir=input_file,
        out_dir=Path(output_file) / 'mod' / lineage,
        file_names=files_to_link,
        overwrite=True,
    )