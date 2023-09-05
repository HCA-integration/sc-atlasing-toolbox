import glob
import yaml
from pprint import pprint
from anndata.experimental import read_elem
import zarr

config = yaml.safe_load(open('test/config.yaml', 'r'))

single_outputs = glob.glob('test/out/preprocessing/*/preprocessed.zarr')
pprint(single_outputs)

for file in single_outputs:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    preprocessing_config = config['DATASETS'][dataset]['preprocessing']
    print(dataset)
    pprint(preprocessing_config)

    print(f'Read {file}...')
    z = zarr.open(file)

    if 'counts' in preprocessing_config['assemble']:
        # assert 'counts' in adata.layers
        X = read_elem(z['X'])
        raw = read_elem(z['raw'])
        assert not isinstance(raw, type(None))
        if X.shape[1] == raw.n_vars:
            assert not all(raw.X.data == X.data)

    if 'normcounts' in preprocessing_config['assemble']:
        assert 'normcounts' in z['layers']

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
