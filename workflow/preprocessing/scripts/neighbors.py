"""
Build kNN graph on embedding
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')

from utils.io import read_anndata, link_zarr
from utils.processing import assert_neighbors
from utils.misc import ensure_dense


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params.get('args', {})
extra_uns = snakemake.params.get('extra_uns', {})

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm', uns='uns', var='var')

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write(output_file)
    exit(0)

# set defaults
if params == False:
    # use existing kNN graph
    adata.obsp = read_anndata(input_file, obs='obs', obsp='obsp').obsp
    assert_neighbors(adata, check_params=False)
    adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
    adata.uns['neighbors']['params'] |= dict(
        connectivities_key='connectivities',
        distances_key='distances',
        use_rep='X',
    )
    assert_neighbors(adata)
else:
    if 'use_rep' not in params:
        if 'X_pca' in adata.obsm:
            params['use_rep'] = 'X_pca'
            params['n_pcs'] = adata.obsm['X_pca'].shape[1]
        else:
            params['use_rep'] = 'X'
    if params['use_rep'] == 'X':
        adata.X = read_anndata(input_file, X='X').X
        ensure_dense(adata)
    logging.info(str(params))

    # compute kNN graph
    sc.pp.neighbors(adata, **params)
    assert_neighbors(adata)

# update .uns
adata.uns |= extra_uns

logging.info(f'Write to {output_file}...')
del adata.X
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obsp', 'uns']]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_keep,
        overwrite=True,
)