import glob
import warnings
from pprint import pformat
import yaml
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


config = yaml.safe_load(open('test/config.yaml', 'r'))

# direct integration outputs
outputs = glob.glob('test/out/marker_genes/dataset~*/file_id~*.zarr')
if len(outputs) == 0:
    warnings.warn('No integration outputs found')

for file in outputs:
    dataset = [x for x in config['DATASETS'] if f'dataset~{x}/' in file][0]

    logging.info(f'dataset: {dataset}')
    logging.info(f'Checking {file}...')
    z = zarr.open(file)

    dataset_config = config['DATASETS'][dataset]['marker_genes']
    logging.info(pformat(dataset_config))
    
    uns = read_elem(z["uns"])
    for group in dataset_config['groups']:
        key = f'marker_genes_group={group}'
        assert key in uns, f'"{key}" not found in {file}\nuns: {uns.keys()}'