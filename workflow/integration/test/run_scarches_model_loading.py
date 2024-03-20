import glob
import warnings
from pprint import pformat
from itertools import chain
import numpy as np
from anndata.experimental import read_elem
from anndata import AnnData
import zarr
import logging
logging.basicConfig(level=logging.INFO)


def get_files_from_pattern(patterns, suffix):
    list_of_lists = [glob.glob(f'{pattern}/{suffix}') for pattern in patterns]
    return list(chain(*list_of_lists))


patterns = [
    f'test/out/integration/dataset~*/file_id~*/batch~*/{pattern}'
    for pattern in ['method~scpoli*', 'method~scgen*']
]
adata_files = get_files_from_pattern(patterns, 'adata.zarr')
model_files = get_files_from_pattern(patterns, 'model')

if len(adata_files) == 0:
    warnings.warn('No integration outputs found')
    exit()

import torch

for adata_file, model_file in zip(adata_files, model_files):
    logging.info(f'Checking {adata_file}...')
    with zarr.open(adata_file) as z:
        adata = AnnData(
            X=read_elem(z["X"]),
            obs=read_elem(z["obs"]),
            var=read_elem(z["var"]),
            uns=read_elem(z["uns"])
        )

    try:
        # check that model history is present
        assert 'model_history' in adata.uns['integration']
        logging.info(pformat(adata.uns['integration']['model_history']))
        
        # check that torch loading works
        model = torch.load(f'{model_file}/model_params.pt')
        logging.info(pformat(list(model.keys())))
    
    except AssertionError as e:
        logging.info('AssertionError for:', adata_file)
        logging.info(pformat(adata.uns))
        raise e
