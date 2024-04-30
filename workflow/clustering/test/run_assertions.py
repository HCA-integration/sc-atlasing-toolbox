import glob
from pprint import pformat
import yaml
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


config = yaml.safe_load(open('test/config.yaml', 'r'))

# direct integration outputss
outputs = glob.glob('test/out/clustering/dataset~*/file_id~*.zarr')
assert len(outputs) > 0, "No outputs found"

for file in outputs:
    dataset = [x for x in config['DATASETS'] if f'dataset~{x}/' in file][0]

    logging.info(f'dataset: {dataset}')
    logging.info(f'Checking {file} for...')
    z = zarr.open(file)

    clustering_config = config['DATASETS'][dataset]['clustering']
    logging.info(pformat(clustering_config))
    
    for res in clustering_config['resolutions']:
        algorithm = clustering_config.get('algorithm', 'leiden')
        cluster_key = f'{algorithm}_{res}'
        assert cluster_key in z["obs"], f"Resolution {cluster_key} not found in {file}"
    