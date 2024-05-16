import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np

from utils.io import read_anndata, get_file_reader


def check_and_set_neighbors_key(adata, neighbors_key):
    neighbors = adata.uns.get(neighbors_key, {})
    neighbors |= {
        'connectivities_key': neighbors.get('connectivities_key', 'connectivities'),
        'distances_key': neighbors.get('distances_key', 'distances'),
        'params': neighbors.get('params', {'use_rep': 'X_pca', 'method': None}),
    }
    adata.uns[neighbors_key] = neighbors
    
    adata.obsp['connectivities'] = adata.obsp[neighbors['connectivities_key']]
    adata.obsp['distances'] = adata.obsp[neighbors['distances_key']]


input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)
neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
cluster_alg = snakemake.params.get('algorithm', 'louvain')
overwrite = snakemake.params.get('overwrite', False)
USE_GPU = False

cluster_key = f'{cluster_alg}_{resolution}'

read_func, _ = get_file_reader(input_file)
with read_func(input_file, 'r') as f:
    cluster_key_exists = cluster_key in f['obs'].keys()

if cluster_key_exists and not overwrite:
    logging.info(f'Read anndata file {input_file} and skip clustering...')
    adata = read_anndata(input_file, obs='obs')
else:
    try:
        import subprocess
        assert subprocess.run('nvidia-smi', shell=True).returncode == 0
        from _rsc_clustering import leiden, louvain
        USE_GPU = True
    except AssertionError:
        logging.info('No GPU found...')
        from scanpy.tools import leiden, louvain

    cluster_alg_map = {
        'louvain': louvain,
        'leiden': leiden,
    }

    logging.info(f'Read anndata file {input_file}...')
    adata = read_anndata(input_file, obs='obs', uns='uns', obsp='obsp')
    # following observations from https://github.com/rapidsai/cugraph/issues/4072#issuecomment-2074822898
    adata.obsp['connectivities'] = adata.obsp['connectivities'].astype('float64')

    # select neighbors
    logging.info(f'Select neighbors for "{neighbors_key}"...')
    check_and_set_neighbors_key(adata, neighbors_key)
    
    logging.info(f'{cluster_alg} clustering with resolution {resolution}...')
    cluster_func = cluster_alg_map.get(cluster_alg, KeyError(f'Unknown clustering algorithm: {cluster_alg}'))
    cluster_func(adata, resolution=resolution, key_added=cluster_key)
    
    if USE_GPU and adata.obs[cluster_key].nunique() > 200:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        import scanpy as sc
        
        alt_cluster_alg_map = {
            'louvain': sc.tl.louvain,
            'leiden': sc.tl.leiden,
        }
        logging.info(f'Cluster {cluster_key} has more than 100 unique values. Falling back to scanpy implementation...')
        cluster_func = alt_cluster_alg_map[cluster_alg]
        cluster_func(adata, resolution=resolution, key_added=cluster_key)

logging.info(f'Write {cluster_key} to {output_file}...')
adata.obs[[cluster_key]].to_parquet(output_file)