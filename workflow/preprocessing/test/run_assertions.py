import warnings
import glob
import yaml
from pprint import pformat
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)

from utils.misc import check_sparse_equal


# config = yaml.safe_load(open('test/config_no_gpu.yaml', 'r'))
# single_outputs = glob.glob('test/out/preprocessing/dataset~no_gpu/file_id~*/preprocessed.zarr')
config = yaml.safe_load(open('test/config.yaml', 'r'))
single_outputs = glob.glob('test/out/preprocessing/dataset~*/file_id~*/preprocessed.zarr')
assert len(single_outputs) > 0, 'No output files found'

for file in single_outputs:
    matching_datasets = [x for x in config['DATASETS'] if x in file]
    if not matching_datasets or matching_datasets[0] == 'empty':
        logging.info(f'No matching dataset found for {file}, skipping...')
        continue
    
    dataset = matching_datasets[0]
    logging.info(f'Check dataset "{dataset}", file {file}')

    preprocessing_config = config['DATASETS'][dataset]['preprocessing']
    logging.info(pformat(preprocessing_config))
    
    z = zarr.open(file)
    try:
        X = read_elem(z['X'])
    except KeyError:
        logging.info(f'No X in {file}, skipping...')
        continue

    if 'normcounts' in preprocessing_config['assemble']:
        assert 'raw' in z, list(z)
        raw = read_elem(z['raw'])
        counts = read_elem(z['layers/counts'])
        assert counts.min() == 0, f'counts.min() = {counts.min()}'
        assert counts.max() > 100, f'counts.max() = {counts.max()}'
        assert not isinstance(raw, type(None))
        if counts.shape[1] == raw.n_vars:
            raw_range = (raw.X.min(), raw.X.max())
            counts_range = (counts.min(), counts.max())
            message = f"matrices aren't equal\nlayers/counts range: {counts_range} vs. raw.X range: {raw_range}"
            assert check_sparse_equal(counts, raw.X), message
            assert any(X.data != raw.X.data)
        
        # check if values are log-transformed
        assert X.min() == 0, f'X.min() = {X.min()}'
        assert X.max() < 20, f'X.max() = {X.max()}'
        
        assert 'normcounts' in z['layers']
        normcounts = read_elem(z['layers/normcounts'])
        assert check_sparse_equal(X, normcounts)

    if 'highly_variable_genes' in preprocessing_config['assemble']:
        var = read_elem(z['var'])
        assert 'highly_variable' in var.columns
        assert 'means' in var.columns
        assert 'dispersions' in var.columns
        assert 'dispersions_norm' in var.columns
        if preprocessing_config['highly_variable_genes'] is False:
            assert var['highly_variable'].sum() == var.shape[0]

    if 'pca' in preprocessing_config['assemble']:
        assert 'X_pca' in z['obsm']

    if 'neighbors' in preprocessing_config['assemble']:
        assert 'neighbors' in z['uns']
        assert 'distances' in z['obsp']
        assert 'connectivities' in z['obsp']
