"""
Assemble multiple preprocessing outputs into single file
"""
from pathlib import Path
from pprint import pformat
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse

from utils.io import read_anndata, link_zarr


def assemble_h5ad(file, file_type, adata):

    def deep_update(d, u):
        """Update metadata from uns
        adapted from:
        https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth#3233356
        """
        for k, v in u.items():
            d[k] = deep_update(d.get(k, {}), v) if isinstance(v, dict) else v
        return d

    logging.info(f'Read {file}...')
    adata_pp = read_anndata(file)
    adata.uns = deep_update(adata.uns, adata_pp.uns)

    if file_type == 'counts':
        logging.debug('add raw counts')
        # adata.layers['counts'] = sparse.csr_matrix(adata_pp[:, adata.var_names].X)
        adata.raw = adata_pp
    elif file_type == 'normalize':
        logging.debug('add normalised counts')
        adata.layers['normcounts'] = sparse.csr_matrix(adata_pp[:, adata.var_names].X)
        adata.X = adata.layers['normcounts']
    elif file_type == 'highly_variable_genes':
        logging.debug('add highly variable gene info')
        hvg_column_map = {
            'highly_variable': False,
            'means': 0,
            'dispersions': 0,
            'dispersions_norm': 0,
            'highly_variable_nbatches': 0,
            'highly_variable_intersection': False,
        }
        hvg_columns = [
            column
            for column in hvg_column_map
            if column in adata_pp.var.columns
        ]
        adata.var[hvg_columns] = adata_pp.var[hvg_columns]
        adata.var = adata.var.fillna(hvg_column_map)
        adata = adata[:, adata_pp.var_names]
    elif file_type == 'pca':
        logging.debug('add PCA')
        adata.obsm['X_pca'] = adata_pp.obsm['X_pca']
        adata.uns['pca'] = adata_pp.obsm['pca']
        adata.uns['varm'] = adata_pp.obsm['varm']
    elif file_type == 'neighbors':
        logging.debug('add neighbors')
        adata.uns['neighbors'] = adata_pp.uns['neighbors']
        adata.obsp['distances'] = adata_pp.obsp['distances']
        adata.obsp['connectivities'] = adata_pp.obsp['connectivities']
    elif file_type == 'umap':
        logging.debug('add UMAP')
        adata.obsm['X_umap'] = adata_pp.obsm['X_umap']
    else:
        ValueError(f'Unknown file type {file_type}')
    return adata


def assemble_zarr(file, file_type, files_to_link):
    file = Path(file)

    if file_type == 'counts':
        logging.debug('add raw counts')
        files_to_link.append((file / 'X', output_file / 'raw' / 'X'))
        files_to_link.append((file / 'obs', output_file / 'obs'))
        files_to_link.append((file / 'var', output_file / 'raw' / 'var'))
    elif file_type == 'normalize':
        logging.debug('add normalised counts')
        files_to_link.append((file / 'X', output_file / 'layers' / 'normcounts'))
        files_to_link.append((file / 'X', output_file / 'X'))
        files_to_link.append((file / 'uns' / 'log1p', output_file / 'uns' / 'log1p'))
        files_to_link.append((file / 'uns' / 'preprocessing' / 'log-transformed', output_file / 'uns' / 'preprocessing' / 'log-transformed'))
        files_to_link.append((file / 'uns' / 'preprocessing' / 'normalization', output_file / 'uns' / 'preprocessing' / 'normalization'))
    elif file_type == 'highly_variable_genes':
        logging.debug('add highly variable gene info')
        files_to_link.append((file / 'X', output_file / 'layers' / 'normcounts'))
        files_to_link.append((file / 'X', output_file / 'X'))
        files_to_link.append((file / 'var', output_file / 'var'))
        files_to_link.append((file / 'uns' / 'preprocessing' / 'highly_variable_genes', output_file / 'uns' / 'preprocessing' / 'highly_variable_genes'))
    elif file_type == 'pca':
        logging.debug('add PCA')
        files_to_link.append((file / 'obsm' / 'X_pca', output_file / 'obsm' / 'X_pca'))
        files_to_link.append((file / 'uns' / 'pca', output_file / 'uns' / 'pca'))
        files_to_link.append((file / 'varm' / 'PCs', output_file / 'varm' / 'PCs'))
    elif file_type == 'neighbors':
        logging.debug('add neighbors')
        files_to_link.append((file / 'obsp' / 'connectivities', output_file / 'obsp' / 'connectivities'))
        files_to_link.append((file / 'obsp' / 'distances', output_file / 'obsp' / 'distances'))
        files_to_link.append((file / 'uns' / 'neighbors', output_file / 'uns' / 'neighbors'))
    elif file_type == 'umap':
        logging.debug('add UMAP')
        files_to_link.append((file / 'obsm' / 'X_umap', output_file / 'obsm' / 'X_umap'))
    else:
        ValueError(f'Unknown file type {file_type}')
    if output_file / 'obs' not in files_to_link:
        files_to_link.append((file / 'obs', output_file / 'obs'))
    if output_file / 'var' not in files_to_link:
        files_to_link.append((file / 'var', output_file / 'raw' / 'var'))
    return files_to_link


output_file = Path(snakemake.output[0])

adata = None
files_to_link = []

for file_type, file in snakemake.input.items():
    if adata is None: # read first file
        logging.info(f'Read first file {file}...')
        adata = read_anndata(file, X='X', obs='obs', var='var')
        if adata.X is None:
            adata.X = np.zeros((adata.n_obs, adata.n_vars))

    if file.endswith('.h5ad'):
        adata = assemble_h5ad(file, file_type, adata)
    elif file.endswith('.zarr'):
        if file_type == 'counts':
            adata.raw = adata
            # adata.layers['counts'] = read_anndata(file).X[:, adata.var_names]
        files_to_link = assemble_zarr(file, file_type, files_to_link)
    else:
        ValueError(f'Unknown file type {file}')

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)

for input_dir, output_dir in files_to_link:
    if not output_dir.exists():
        logging.info(f'Directory {output_dir} does not yet exist, create it...')
        output_dir.mkdir(parents=True)
    link_zarr(
        in_dir=input_dir,
        out_dir=output_dir,
        overwrite=True,
    )