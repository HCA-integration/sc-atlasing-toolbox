import mudata as mu
import logging
logging.basicConfig(level=logging.INFO)

from metrics.utils import anndata_to_mudata, compute_neighbors, get_from_adata
from utils.io import read_anndata_or_mudata


input_adata = snakemake.input[0]
output_file = snakemake.output[0]
lineage_key = snakemake.wildcards.lineage_key

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

logging.info(mudata.__str__())
mudata.write_zarr(output_file)