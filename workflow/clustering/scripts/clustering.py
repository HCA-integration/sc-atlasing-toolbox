import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np

from utils.io import read_anndata, write_zarr_linked, get_file_reader
from utils.misc import dask_compute


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


def apply_clustering(
    adata,
    resolution: float,
    cpu_kwargs: dict = None,
    n_cell_cpu: int = 300_000,
    max_clusters: bool = None,
    recompute_neighbors: bool = False,
    neighbors_args: dict = {},
    **kwargs,
):
    """
    :param adata: anndata object
    :param cpu_kwargs: clustering parameters for CPU implementation
    :param n_cell_cpu: number of cells to use CPU implementation
    :param max_clusters: maximum number of clusters to determine if number of clusters are correct
    """
    import scanpy as sc
    USE_GPU = False
    try:
        import subprocess
        assert subprocess.run('nvidia-smi', shell=True, stdout=subprocess.DEVNULL).returncode == 0
        from rapids_singlecell.tl import leiden, louvain
        from rapids_singlecell.pp import neighbors
        USE_GPU = True
    except Exception as e:
        logging.info(f'Error importing rapids found...\n{e}')
        from scanpy.tools import leiden, louvain
        from scanpy.preprocessing import neighbors
    
    algorithm_map = {
        'louvain': louvain,
        'leiden': leiden,
    }
    alt_algorithm_map = {
        'louvain': sc.tl.louvain,
        'leiden': sc.tl.leiden,
    }
    
    kwargs |= dict(resolution=resolution)
    
    if not cpu_kwargs:
        cpu_kwargs = dict()
    
    if not USE_GPU:
        kwargs |= cpu_kwargs
    
    if recompute_neighbors:
        logging.info(f'Compute neighbors for {neighbors_args}...')
        neighbors(adata, **neighbors_args)
    
    if USE_GPU:
        # following observations from https://github.com/rapidsai/cugraph/issues/4072#issuecomment-2074822898
        adata.obsp['connectivities'] = adata.obsp['connectivities'].astype('float64')

    if adata.n_obs < n_cell_cpu:
        # switch to CPU implementation for smaller numbers of cells
        algorithm_map = alt_algorithm_map
        
    logging.info(f'{algorithm} clustering with {kwargs} for {adata.n_obs} cells...')
    cluster_func = algorithm_map.get(algorithm, KeyError(f'Unknown clustering algorithm: {algorithm}'))
    cluster_func(adata, **kwargs)
    
    if not max_clusters:
        max_clusters = max(1, int(50 * resolution))
    n_clusters = adata.obs[cluster_key].nunique()
    
    if USE_GPU and n_clusters > max_clusters:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        logging.info(
            f'Cluster {cluster_key} has {n_clusters} custers, which is more than {max_clusters}. '
            'Falling back to scanpy implementation...'
        )
        cluster_func = alt_algorithm_map[algorithm]
        kwargs |= cpu_kwargs
        cluster_func(adata, **kwargs)
    return adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)
algorithm = snakemake.wildcards.algorithm
level = int(snakemake.wildcards.level)
threads = snakemake.threads
overwrite = snakemake.params.get('overwrite', False)

# set parameters for clustering
cluster_key = f'{algorithm}_{resolution}_{level}'
kwargs = dict(
    resolution=resolution,
    key_added=cluster_key,
) | snakemake.params.get('clustering_args', {})
cpu_kwargs = dict(flavor='igraph')
if algorithm == 'leiden':
    cpu_kwargs |= dict(n_iterations=2)

read_func, _ = get_file_reader(input_file)
with read_func(input_file, 'r') as f:
    cluster_key_exists = cluster_key in f['obs'].keys()

if cluster_key_exists and not overwrite:
    logging.info(f'Read anndata file {input_file} and skip clustering...')
    adata = read_anndata(input_file, obs='obs')
else:
    read_kwargs = dict(obs='obs', obsm='obsm')
    
    if level <= 1:
        logging.info(f'Read anndata file {input_file}...')
        if not input_file.endswith('.h5ad'):
            read_kwargs.pop('obsm', None) # not needed, since neighobrs are already computed
        adata = read_anndata(input_file, uns='uns', obsp='obsp', **read_kwargs)
        
        # select neighbors
        neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
        logging.info(f'Select neighbors for "{neighbors_key}"...')
        check_and_set_neighbors_key(adata, neighbors_key)
                
        adata = apply_clustering(
            adata,
            cpu_kwargs=cpu_kwargs,
            recompute_neighbors=False,
            **kwargs,
        )
    else:
        from concurrent.futures import ProcessPoolExecutor
        from concurrent.futures import as_completed
        
        logging.info(f'Read anndata file {input_file}...')
        adata = read_anndata(input_file, dask=True, backed=True, **read_kwargs)
        
        def cluster_subset(
            adata,
            resolution,
            key_added,
            prev_cluster_key,
            prev_cluster_value,
            neighbors_args, # TODO: custom parameters
        ):
            logging.info(f'Subsetting to {prev_cluster_key}={prev_cluster_value}...')
            sub_adata = dask_compute(adata[adata.obs[prev_cluster_key] == prev_cluster_value].copy())
            
            if sub_adata.n_obs < 2 * neighbors_args.get('n_neighbors', 15):
                prev_cluster_clean = sub_adata.obs[prev_cluster_key].astype(str).apply(lambda x: x.split('_')[-1])
                return sub_adata.obs[prev_cluster_key].str.cat(prev_cluster_clean, sep='_')
            
            sub_adata = apply_clustering(
                sub_adata,
                resolution=resolution,
                key_added=key_added,
                cpu_kwargs=cpu_kwargs,
                neighbors_args=neighbors_args,
                recompute_neighbors=True,
            )
            return sub_adata.obs[[prev_cluster_key, key_added]].agg('_'.join, axis=1)
        
        
        prev_cluster_key = f'{algorithm}_{resolution}_{level-1}'
        neighbors_args = dict(use_rep='X_pca') | snakemake.params.get('neighbors_args', {})
        
        logging.info(f'Cluster {cluster_key} with {threads} threads...')
        with ProcessPoolExecutor(threads) as executor:
            futures = [
                executor.submit(
                    cluster_subset,
                    adata,
                    prev_cluster_key=prev_cluster_key,
                    prev_cluster_value=prev_cluster,
                    neighbors_args=neighbors_args,
                    **kwargs
                )
                for prev_cluster in adata.obs[prev_cluster_key].unique()
            ]

            for future in as_completed(futures):
                clusters = future.result()
                adata.obs.loc[clusters.index, cluster_key] = clusters


logging.info(f'Write {cluster_key} to {output_file}...')
adata.obs = adata.obs[[cluster_key]]
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)