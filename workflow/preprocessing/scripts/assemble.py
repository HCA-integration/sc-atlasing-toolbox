"""
Assemble multiple preprocessing outputs into single file
"""
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse

from utils.io import read_anndata, link_zarr


def link_file(input_dir, output_dir, files_to_keep=None):
    if not output_dir.exists():
        logging.info(f'Create directory {output_dir}')
        output_dir.mkdir(parents=True)
    if files_to_keep is None:
        files_to_keep = []
    input_files = [f.name for f in Path(input_dir).iterdir()]
    link_zarr(
        in_dir=input_dir,
        out_dir=output_dir,
        file_names=[f for f in input_files if f not in files_to_keep],
        overwrite=True,
    )


def deep_update(d, u):
    """Update metadata from uns
    adapted from:
    https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth#3233356
    """
    for k, v in u.items():
        d[k] = deep_update(d.get(k, {}), v) if isinstance(v, dict) else v
    return d


output_file = Path(snakemake.output[0])

adata = None

# TODO: link to other zarr directories

for file_type, file in snakemake.input.items():
    logging.info(f'Read {file}...')
    adata_pp = read_anndata(file)

    if adata is None:
        adata = adata_pp
        logging.info(f'Write to {output_file}...')
        adata.write_zarr(output_file)
        
    if file.endswith('.h5ad'):
        if not output_file.exists():
            adata.write_zarr(output_file)
        file = output_file
    file = Path(file)

    if file_type == 'counts':
        logging.debug('add raw counts')
        adata.layers['counts'] = sparse.csr_matrix(adata_pp[:, adata.var_names].X)
        link_file(file / 'X', output_file / 'layers' / 'counts')

    elif file_type == 'normalize':
        logging.debug('add normalised counts')
        adata.layers['normcounts'] = sparse.csr_matrix(adata_pp[:, adata.var_names].X)
        adata.X = adata.layers['normcounts']
        link_file(file / 'X', output_file / 'layers' / 'normcounts')
        link_file(file / 'X', output_file / 'X')
        link_file(file / 'uns' / 'log1p', output_file / 'uns' / 'log1p')
        link_file(file / 'uns' / 'preprocessing' / 'log-transformed', output_file / 'uns' / 'preprocessing' / 'log-transformed')
        link_file(file / 'uns' / 'preprocessing' / 'normalization', output_file / 'uns' / 'preprocessing' / 'normalization')

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
        link_file(file / 'var', output_file / 'var')
        link_file(file / 'uns' / 'preprocessing' / 'highly_variable_genes', output_file / 'uns' / 'preprocessing' / 'highly_variable_genes')

    elif file_type == 'pca':
        logging.debug('add PCA')
        adata.obsm['X_pca'] = adata_pp.obsm['X_pca']
        link_file(file / 'obsm' / 'X_pca', output_file / 'obsm' / 'X_pca')
        link_file(file / 'uns' / 'pca', output_file / 'uns' / 'pca')
        link_file(file / 'varm' / 'PCs', output_file / 'varm' / 'PCs')

    elif file_type == 'neighbors':
        logging.debug('add neighbors')
        adata.uns['neighbors'] = adata_pp.uns['neighbors']
        adata.obsp['distances'] = adata_pp.obsp['distances']
        adata.obsp['connectivities'] = adata_pp.obsp['connectivities']
        link_file(file / 'obsp' / 'connectivities', output_file / 'obsp' / 'connectivities')
        link_file(file / 'obsp' / 'distances', output_file / 'obsp' / 'distances')
        link_file(file / 'uns' / 'neighbors', output_file / 'uns' / 'neighbors')

    elif file_type == 'umap':
        logging.debug('add UMAP')
        adata.obsm['X_umap'] = adata_pp.obsm['X_umap']
        link_file(file / 'obsm' / 'X_umap', output_file / 'obsm' / 'X_umap')

    else:
        ValueError(f'Unknown file type {file_type}')

    adata.uns = deep_update(adata.uns, adata_pp.uns)

logging.info(adata.__repr__())
logging.info(list(output_file.iterdir()))

# logging.info('Preprocessing metadata in adata.uns:')
# logging.info(pformat(adata.uns['preprocessing']))
