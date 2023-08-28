import anndata
import glob
import yaml
from pprint import pprint

config = yaml.safe_load(open('test/config.yaml', 'r'))

single_outputs = glob.glob('test/out/preprocessing/*/preprocessed.zarr')
print(single_outputs)

for file in single_outputs:
    adata = anndata.read_zarr(file)
    dataset = [x for x in config['DATASETS'] if x in file][0]
    preprocessing_config = config['DATASETS'][dataset]['preprocessing']

    print(dataset)
    print(file)
    pprint(preprocessing_config)

    if 'counts' in preprocessing_config['assemble']:
        # assert 'counts' in adata.layers
        assert not isinstance(adata.raw, type(None))
        if adata.n_vars == adata.raw.n_vars:
            assert not all(adata.raw.X.data == adata.X.data)

    if 'normcounts' in preprocessing_config['assemble']:
        assert 'normcounts' in adata.layers

    if 'highly_variable_genes' in preprocessing_config['assemble']:
        assert 'highly_variable' in adata.var.columns
        assert 'means' in adata.var.columns
        assert 'dispersions' in adata.var.columns
        assert 'dispersions_norm' in adata.var.columns
        if preprocessing_config['highly_variable_genes'] is False:
            assert adata.var['highly_variable'].sum() == adata.n_vars

    if 'pca' in preprocessing_config['assemble']:
        assert 'X_pca' in adata.obsm

    if 'neighbors' in preprocessing_config['assemble']:
        assert 'neighbors' in adata.uns
        assert 'distances' in adata.obsp
        assert 'connectivities' in adata.obsp
