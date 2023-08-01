import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from methods.utils import read_anndata
from metrics.utils import compute_neighbors, get_from_adata


input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)
meta = get_from_adata(adata)

logging.info('Compute neighbors...')
for output_type in meta['output_types']:
    logging.info(f'Computing neighbors for output type {output_type}...')
    # compute neighbors by output type
    compute_neighbors(adata, output_type)
    adata.obsp[f'connectivities_{output_type}'] = adata.obsp['connectivities']
    adata.obsp[f'distances_{output_type}'] = adata.obsp['distances']

    logging.info(f'Leiden clustering with resolution {resolution}...')
    cluster_key = f'leiden_{resolution}_{output_type}'
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

logging.info('Write file...')
adata.obs[cluster_key].to_csv(output_dir, sep='\t', index=True)
