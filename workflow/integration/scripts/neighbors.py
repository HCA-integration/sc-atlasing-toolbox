import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from metrics.utils import compute_neighbors, get_from_adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)
meta = get_from_adata(adata)

for output_type in meta['output_types']:
    logging.info(f'Computing neighbors for output type {output_type}...')
    # compute neighbors by output type
    adata_knn = compute_neighbors(adata, output_type)
    adata.obsp[f'connectivities_{output_type}'] = adata_knn.obsp['connectivities']
    adata.obsp[f'distances_{output_type}'] = adata_knn.obsp['distances']

logging.info('Write file...')
del adata.X
del adata.layers
del adata.obsm
adata.write_zarr(output_file)
