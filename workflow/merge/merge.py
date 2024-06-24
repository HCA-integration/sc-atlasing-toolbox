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
from utils.misc import apply_layers, dask_compute


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

# subset to non-empty datasets
files = [file for file, adata in zip(files, adatas) if adata.n_obs > 0]
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
    
    # concatenate
    adata = sc.concat(adatas, join=merge_strategy)
    
elif backed:
    logging.info('Read all files in backed mode...')
    adatas = [read_adata(file, backed, dask) for file in files]
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
        logging.info(f'{file}:\n{_adata}')
        
        # merge adata
        adata = sc.concat([adata, _adata], join=merge_strategy)
        logging.info(f'after merge:\n{adata}')

        del _adata
        gc.collect()


# merge lost annotations
logging.info('Add gene info...')
for _ad in adatas:
    # get intersection of var_names
    var_names = list(set(adata.var_names).intersection(set(_ad.var_names)))
    
    # conver categorical columns to str
    categorical_columns = _ad.var.select_dtypes(include='category').columns
    _ad.var[categorical_columns] = _ad.var[categorical_columns].astype(str)
    
    # add new columns to adata.var
    adata.var.loc[var_names, _ad.var.columns] = _ad.var.loc[var_names, :]

# fix dtypes
adata.var = adata.var.infer_objects()
logging.info(adata.var)

# set new indices
adata.obs_names = dataset + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# add uns data
adata.uns['dataset'] = dataset
logging.info(adata.__str__())

logging.info('Compute matrix...')
adata = dask_compute(adata)

logging.info(f'Write to {out_file}...')
adata.write_zarr(out_file)