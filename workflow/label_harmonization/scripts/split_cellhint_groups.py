import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path

from utils.io import read_anndata


input_file = snakemake.input[0]
output_dir = snakemake.output[0]
Path(output_dir).mkdir(parents=True)

adata = read_anndata(input_file, obs='obs')
groups = adata.obs['group'].unique().dropna().tolist()

logging.info(f'Write groups to {output_dir}...')
for group in groups+['all']:
    group_file = f'{output_dir}/{group}.yaml'
    logging.info(group_file)
    open(group_file, 'w').close()