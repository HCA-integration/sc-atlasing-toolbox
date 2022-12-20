from pathlib import Path
import gc

import pandas as pd
import scanpy as sc
import anndata

from utils import CELLxGENE_VARS


def read_adata(file):
    ad = anndata.read_zarr(file)
    ad.var = ad.var[CELLxGENE_VARS]
    # remove data
    del ad.uns
    del ad.layers
    del ad.raw
    del ad.obsm
    return ad


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.zarr

if len(files) == 1:
    in_file = Path(files[0])
    out_dir = Path(out_file)
    if not out_dir.exists():
        out_dir.mkdir()
    for f in in_file.iterdir():
        if f.name == '.snakemake_timestamp':
            continue  # skip snakemake timestamp
        new_file = out_dir / f.name
        new_file.symlink_to(f.resolve())
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
    adata.write_zarr(out_file)
