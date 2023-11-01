import logging
logging.basicConfig(level=logging.INFO)
import gc

import pandas as pd
import scanpy as sc
import anndata

from utils import SCHEMAS, get_union
from utils_pipeline.io import link_zarr


def read_adata(file, keep_columns):
    ad = anndata.read_zarr(file)
    if not keep_columns:
        logging.info('Keep only mandatory columns...')
        ad.obs = ad.obs[get_union(SCHEMAS["CELLxGENE_OBS"], SCHEMAS["EXTRA_COLUMNS"])]
    ad.var = ad.var[SCHEMAS['CELLxGENE_VARS']]
    # remove data
    del ad.layers
    del ad.raw
    del ad.obsm
    return ad


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.zarr
merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
else:
    logging.info(f'Read first file {files[0]}...')
    adata = read_adata(files[0], keep_all_columns)
    logging.info(adata.__str__())

    uns_per_dataset = {
        adata.uns['meta']['dataset']: adata.uns['meta']
    }

    for file in files[1:]:
        logging.info(f'Read {file}...')
        _adata = read_adata(file, keep_all_columns)
        logging.info(_adata.__str__())

        if _adata.n_obs == 0:
            logging.info('Empty adata, skip concatenation...')
            continue

        # collect metadata
        uns_per_dataset[_adata.uns['meta']['dataset']] = _adata.uns['meta']

        logging.info('Concatenate...')
        # merge genes
        var_map = pd.merge(
            adata.var,
            _adata.var,
            how=merge_strategy,
            on=['feature_id'] + SCHEMAS['CELLxGENE_VARS']
        )
        var_map = var_map[~var_map.index.duplicated()]

        # merge adata
        adata = sc.concat([adata, _adata], join='outer')
        adata = adata[:, var_map.index]
        adata.var = var_map.loc[adata.var_names]

        del _adata
        gc.collect()

    organ = adata.obs['organ'].unique()
    assert len(organ) == 1
    adata.uns['dataset'] = dataset
    adata.uns['organ'] = organ
    adata.uns['meta'] = {
        'dataset': dataset,
        'organ': organ,
        'per_dataset': uns_per_dataset,
    }
    print(adata)
    print(adata.var)

    assert 'feature_name' in adata.var

    logging.info('Write...')
    adata.write_zarr(out_file)
