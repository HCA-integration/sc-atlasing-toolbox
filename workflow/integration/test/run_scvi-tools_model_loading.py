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
    for pattern in ['method~sc*vi*']
]
adata_files = get_files_from_pattern(patterns, 'adata.zarr')
model_files = get_files_from_pattern(patterns, 'model')

if len(adata_files) == 0:
    warnings.warn('No integration outputs found')
    exit()

import scvi
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
        model = torch.load(f'{model_file}/model.pt')
        
        # check that scVI loading works
        method = adata.uns['integration']['method']
        if method == 'scvi':
            model = scvi.model.SCVI.load(model_file, adata)
        elif method == 'scanvi':
            model = scvi.model.SCANVI.load(model_file, adata)
        else:
            raise ValueError(f'Unknown model for {adata_file}')
        logging.info(pformat(model))
    
    except AssertionError as e:
        logging.info('AssertionError for:', adata_file)
        logging.info(pformat(adata.uns))
        raise e
