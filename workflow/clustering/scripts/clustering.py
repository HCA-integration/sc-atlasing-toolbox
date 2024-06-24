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

# set parameters for clustering
cluster_key = f'{cluster_alg}_{resolution}'
kwargs = dict(resolution=resolution, key_added=cluster_key)
cpu_kwargs = dict(flavor='igraph')
if cluster_alg == 'leiden':
    cpu_kwargs |= dict(n_iterations=2)

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
    except Exception as e:
        logging.info(f'Error importing rapids found...\n{e}')
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
    
    if not USE_GPU:
        kwargs |= cpu_kwargs
    logging.info(f'{cluster_alg} clustering with {kwargs}...')
    cluster_func = cluster_alg_map.get(cluster_alg, KeyError(f'Unknown clustering algorithm: {cluster_alg}'))
    cluster_func(adata, **kwargs)
    
    max_clusters = max(1, int(50 * resolution))
    n_clusters = adata.obs[cluster_key].nunique()
    if USE_GPU and n_clusters > max_clusters:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        import scanpy as sc
        
        alt_cluster_alg_map = {
            'louvain': sc.tl.louvain,
            'leiden': sc.tl.leiden,
        }
        logging.info(
            f'Cluster {cluster_key} has {n_clusters} custers, which is more than {max_clusters}.'
            'Falling back to scanpy implementation...'
        )
        cluster_func = alt_cluster_alg_map[cluster_alg]
        kwargs |= cpu_kwargs
        cluster_func(adata, **kwargs)

logging.info(f'Write {cluster_key} to {output_file}...')
adata.obs[[cluster_key]].to_parquet(output_file)