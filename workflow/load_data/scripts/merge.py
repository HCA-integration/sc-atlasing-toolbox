from pathlib import Path
import gc

import pandas as pd
import scanpy as sc

from utils import CELLxGENE_VARS


def read_adata(file):
    ad = sc.read(file)
    ad.var = ad.var[CELLxGENE_VARS]
    # remove data
    del ad.uns
    del ad.layers
    del ad.raw
    del ad.obsm
    return ad


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.h5ad

if len(files) == 1:
    Path(out_file).symlink_to(Path(files[0]).resolve())
else:
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
    adata.write(out_file, compression='gzip')
