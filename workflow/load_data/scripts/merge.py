import logging
logging.basicConfig(level=logging.INFO)
import gc

import pandas as pd
import scanpy as sc
from anndata.experimental import AnnCollection
from anndata import AnnData

from utils import SCHEMAS, get_union
from utils_pipeline.io import read_anndata, link_zarr


def read_adata(file, keep_columns, backed=True):
    logging.info(f'Read {file}...')
    ad = read_anndata(file, backed=backed, X='X', obs='obs', var='var', uns='uns')
    if not keep_columns:
        logging.info(f'Keep only mandatory columns for {file}...')
        ad.obs = ad.obs[get_union(SCHEMAS["CELLxGENE_OBS"], SCHEMAS["EXTRA_COLUMNS"])]
    ad.var = ad.var[SCHEMAS['CELLxGENE_VARS']]
    return ad


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.zarr

# merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)
backed = snakemake.params.get('backed', True)

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
    exit(0)

adatas = [read_adata(file, keep_all_columns, backed) for file in files]
print(adatas)
adatas = [adata for adata in adatas if adata.n_obs > 0]

if len(adatas) == 0:
    logging.info('All adatas are empty, skip concatenation...')
    AnnData().write_zarr(out_file)
    exit(0)

dc = AnnCollection(
    adatas,
    join_obs='outer',
    join_obsm=None,
    join_vars='inner',
    indices_strict=not backed,
)
logging.info(dc.__str__())

logging.info('Subset AnnDataCollection, returning a View...')
adata = dc[:].to_adata()
assert adata.X is not None

# set new indices
adata.obs_names = dataset + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# merge lost annotations
var_map = adatas[0].var
uns = adatas[0].uns
uns_per_dataset = {
    uns['meta']['dataset']: uns['meta']
}
for _adata in adatas[1:]:
    # uns
    uns_per_dataset[_adata.uns['meta']['dataset']] = _adata.uns['meta']
    
    # genes
    var_map = pd.merge(
        var_map,
        _adata.var,
        how='inner',
        on=['feature_id'] + SCHEMAS['CELLxGENE_VARS']
    )
    var_map = var_map[~var_map.index.duplicated()]

logging.info('Add gene info...')
adata = adata[:, var_map.index]
adata.var = var_map.loc[adata.var_names]
logging.info(adata.var)
assert 'feature_name' in adata.var

# add uns data
organ = adata.obs['organ'].unique()
assert len(organ) == 1
adata.uns['dataset'] = dataset
adata.uns['organ'] = organ
adata.uns['meta'] = {
    'dataset': dataset,
    'organ': organ,
    'per_dataset': uns_per_dataset,
}
logging.info(adata.__str__())

logging.info(f'Write to {out_file}...')
adata.write_zarr(out_file)


# logging.info(f'Read first file {files[0]}...')
# adata = read_adata(files[0], keep_all_columns)
# logging.info(adata.__str__())

# uns_per_dataset = {
#     adata.uns['meta']['dataset']: adata.uns['meta']
# }

# for file in files[1:]:
#     logging.info(f'Read {file}...')
#     _adata = read_adata(file, keep_all_columns)
#     logging.info(_adata.__str__())

#     if _adata.n_obs == 0:
#         logging.info('Empty adata, skip concatenation...')
#         continue

#     # collect metadata
#     uns_per_dataset[_adata.uns['meta']['dataset']] = _adata.uns['meta']

#     logging.info('Concatenate...')
#     # merge genes
#     var_map = pd.merge(
#         adata.var,
#         _adata.var,
#         how=merge_strategy,
#         on=['feature_id'] + SCHEMAS['CELLxGENE_VARS']
#     )
#     var_map = var_map[~var_map.index.duplicated()]

#     # merge adata
#     adata = sc.concat([adata, _adata], join='outer')
#     adata = adata[:, var_map.index]
#     adata.var = var_map.loc[adata.var_names]

#     del _adata
#     gc.collect()
