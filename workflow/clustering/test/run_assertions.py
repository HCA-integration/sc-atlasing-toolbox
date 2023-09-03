import glob
from pprint import pprint
import yaml
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


config = yaml.safe_load(open('test/config.yaml', 'r'))

# direct integration outputss
outputs = glob.glob('test/out/clustering/resolutions/*.zarr')
pprint(outputs)

for file in outputs:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    clustering_config = config['DATASETS'][dataset]['preprocessing']
    print(dataset)
    pprint(clustering_config)

    logging.info(f'Checking {file}...')
    z = zarr.open(file)
    # obs = read_elem(z["obs"])

    try:
        for res in clustering_config['resolutions']:
            assert res in z["obs"]
    except AssertionError as e:
        print('AssertionError for:', file)
        print(z["obs"])
        raise e
