import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from metrics.utils import get_from_adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)
meta = get_from_adata(adata)

for output_type in meta['output_types']:
    adata.obsp['connectivities'] = adata.obsp[f'connectivities_{output_type}']
    adata.obsp['distances'] = adata.obsp[f'distances_{output_type}']

    logging.info(f'Leiden clustering with resolution {resolution} and output type {output_type}...')
    cluster_key = f'leiden_{resolution}_{output_type}'
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

logging.info('Write file...')
adata.obs[cluster_key].to_csv(output_file, sep='\t', index=True)
