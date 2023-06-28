import scanpy as sc
import glob
import yaml
from pprint import pprint

config = yaml.safe_load(open('test/config.yaml', 'r'))
pprint(config['DATASETS'])

single_outputs = glob.glob('test/out/preprocessing/*/preprocessed.h5ad')
print(single_outputs)

for file in single_outputs:
    adata = sc.read(file)
    dataset = [x for x in config['DATASETS'] if x in file][0]
    preprocessing_config = config['DATASETS'][dataset]['preprocessing']['assemble']
    
    try:
        if 'counts' in preprocessing_config:
            assert 'counts' in adata.layers
        if 'normcounts' in preprocessing_config:
            assert 'normcounts' in adata.layers
        if 'highly_variable_genes' in preprocessing_config:
            assert 'highly_variable' in adata.var.columns
            assert 'means' in adata.var.columns
            assert 'dispersions' in adata.var.columns
            assert 'dispersions_norm' in adata.var.columns
        if 'pca' in preprocessing_config:
            assert 'X_pca' in adata.obsm
        if 'neighbors' in preprocessing_config:
            assert 'neighbors' in adata.uns
            assert 'distances' in adata.obsp
            assert 'connectivities' in adata.obsp
    except AssertionError as e:
        print(dataset)
        pprint(preprocessing_config)
        raise e