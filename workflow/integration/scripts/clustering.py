import sys
from pathlib import Path
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from methods.utils import read_anndata


input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)

logging.info(f'Leiden clustering with resolution {resolution}...')
cluster_key = f'leiden_{resolution}'
sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

logging.info('Write file...')
adata.obs[cluster_key].to_csv(output_dir, sep='\t', index=True)
