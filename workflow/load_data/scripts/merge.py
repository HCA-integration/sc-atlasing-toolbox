import gc

import pandas as pd
import scanpy as sc

from utils import CELLxGENE_OBS, CELLxGENE_VARS, get_union

dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.h5ad


def read_adata(file):
    ad = sc.read(file)
    ad.obs_names = ad.uns['dataset'] + '-' + ad.obs.reset_index().index.astype(str)

    # keep only relevant columns
    extra_columns = ['organ', 'donor', 'sample', 'cell_annotation', 'reference', 'study', 'dataset', 'modalities']
    columns = get_union(CELLxGENE_OBS, extra_columns)

    ad.obs = ad.obs[columns].copy()

    ad.var = ad.var[CELLxGENE_VARS]
    ad.var.index.set_names('feature_id', inplace=True)

    # remove data
    del ad.uns
    del ad.layers
    del ad.raw
    del ad.obsm

    gc.collect()

    return ad


print(f'Read first file {files[0]}...')
adata = read_adata(files[0])
print(adata)

for file in files[1:]:
    print(f'Read {file}...')
    _adata = read_adata(file)
    print(_adata)

    print('Concatenate...')
    # merge genes
    var_map = pd.merge(
        adata.var,
        _adata.var,
        how='outer',
        on=['feature_id'] + CELLxGENE_VARS
    )
    var_map = var_map[~var_map.index.duplicated()]

    # merge adata
    adata = sc.concat([adata, _adata], join='outer')
    adata.var = var_map.loc[adata.var_names]

    del _adata
    gc.collect()

organ = adata.obs['organ'].unique()[0]
adata.uns['dataset'] = dataset
adata.uns['organ'] = organ
adata.uns['meta'] = {'dataset': dataset, 'organ': organ}
print(adata)
print(adata.var)

assert 'feature_name' in adata.var

print('Write...')
adata.write(out_file)
