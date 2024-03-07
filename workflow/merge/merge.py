import logging
logging.basicConfig(level=logging.INFO)
import gc
import faulthandler
faulthandler.enable()
import pandas as pd
import scanpy as sc
from anndata.experimental import AnnCollection
from anndata import AnnData

from utils.io import read_anndata, link_zarr
from utils.misc import apply_layers


def read_adata(file, backed=False, dask=False):
    logging.info(f'Read {file}...')
    adata = read_anndata(
        file,
        backed=backed,
        dask=dask,
        X='X',
        obs='obs',
        var='var',
        uns='uns',
    )
    return adata


dataset = snakemake.wildcards.dataset
files = snakemake.input
out_file = snakemake.output.zarr

merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
    exit(0)

adatas = [read_anndata(file, obs='obs', var='var', uns='uns') for file in files]
adatas = [adata for adata in adatas if adata.n_obs > 0]

if len(adatas) == 0:
    logging.info('All adatas are empty, skip concatenation...')
    AnnData().write_zarr(out_file)
    exit(0)

if dask:
    from dask import array as da
    from dask import config as da_config

    da_config.set(
        **{
            'num_workers': snakemake.threads,
            'array.slicing.split_large_chunks': True
        }
    )
    logging.info('Read all files with dask...')
    logging.info(f'n_threads: {snakemake.threads}')
    adatas = [read_adata(file, backed, dask) for file in files]
    adatas = [adata for adata in adatas if adata.n_obs > 0]
    
    # concatenate
    adata = sc.concat(adatas, join=merge_strategy)
    
elif backed:
    logging.info('Read all files in backed mode...')
    adatas = [read_adata(file, backed, dask) for file in files]
    adatas = [adata for adata in adatas if adata.n_obs > 0]
    dc = AnnCollection(
        adatas,
        join_obs='outer',
        join_obsm=None,
        join_vars=merge_strategy,
        indices_strict=not backed,
    )
    logging.info(dc.__str__())

    logging.info('Subset AnnDataCollection, returning a View...')
    adata = dc[:].to_adata()
    assert adata.X is not None

else:
    logging.info(f'Read first file {files[0]}...')
    adata = read_adata(files[0], backed=backed, dask=dask)
    logging.info(adata.__str__())

    for file in files[1:]:
        logging.info(f'Read {file}...')
        _adata = read_adata(file, backed=backed, dask=dask)
        logging.info(_adata.__str__())
        
        # merge adata
        adata = sc.concat([adata, _adata], join=merge_strategy)

        del _adata
        gc.collect()

# set new indices
adata.obs_names = dataset + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# add uns data
adata.uns['dataset'] = dataset
logging.info(adata.__str__())

logging.info(f'Write to {out_file}...')
adata.write_zarr(out_file)