"""
Assemble multiple preprocessing outputs into single file
"""
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from scipy import sparse

from utils.io import read_anndata

output_file = snakemake.output[0]

adata = None

# TODO: link to other zarr directories

for file_type, file in snakemake.input.items():
    logging.info(f'Read {file}...')
    adata_pp = read_anndata(file)
    
    if adata is None:
        adata = adata_pp
    
    if file_type == 'counts':
        logging.debug('add raw counts')
        adata.layers['counts'] = sparse.csr_matrix(adata_pp[:, adata.var_names].X)

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
        hvg_columns = [column for column in hvg_column_map.keys() if column in adata_pp.var.columns]
        adata.var[hvg_columns] = adata_pp.var[hvg_columns]
        adata.var = adata.var.fillna(hvg_column_map)

    elif file_type == 'pca':
        logging.debug('add PCA')
        adata.obsm['X_pca'] = adata_pp.obsm['X_pca']

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


    # update metadata from uns
    def update(d, u):
        """
        adapted from:
        https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth#3233356
        """
        for k, v in u.items():
            d[k] = update(d.get(k, {}), v) if isinstance(v, dict) else v
        return d


    adata.uns = update(adata.uns, adata_pp.uns)

logging.info(adata.__repr__())
logging.info('Preprocessing metadata in adata.uns:')
logging.info(pformat(adata.uns['preprocessing']))

logging.info(f'Write to {output_file}...')
adata.write_zarr(output_file)
