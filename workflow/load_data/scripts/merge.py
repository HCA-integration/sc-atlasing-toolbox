import gc

import pandas as pd
import scanpy as sc

from utils import CELLxGENE_OBS, CELLxGENE_VARS, get_union

organ = snakemake.wildcards.organ
files = snakemake.input
out_file = snakemake.output.h5ad


def read_adata(file):
    ad = sc.read(file)
    ad.obs_names = ad.uns['meta']['dataset_name'] + '-' + ad.obs.reset_index().index.astype(str)

    # keep only relevant columns
    columns = get_union(CELLxGENE_OBS, ['donor', 'sample', 'cell_annotation', 'reference'])

    ad.obs = ad.obs[columns].copy()
    ad.obs['organ'] = organ
    ad.obs['dataset'] = ad.uns['meta']['dataset_name']
    ad.obs['dataset_id'] = ad.uns['meta']['dataset_id']

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
        on=['feature_id']+CELLxGENE_VARS
    )
    var_map = var_map[~var_map.index.duplicated()]

    # merge adata
    adata = sc.concat([adata, _adata], join='outer')
    adata.var = var_map.loc[adata.var_names]

    del _adata
    gc.collect()

adata.uns['dataset'] = organ
adata.uns['organ'] = organ
print(adata)
print(adata.var)

assert 'feature_name' in adata.var

print('Write...')
adata.write(out_file)
