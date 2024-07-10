import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np

from utils.io import read_anndata, write_zarr_linked, get_file_reader


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
    use_gpu=False,
    cpu_kwargs={},
    max_clusters=None,
    **kwargs,
):
    import scanpy as sc
    
    algorithm_map = {
        'louvain': louvain,
        'leiden': leiden,
    }
    alt_algorithm_map = {
        'louvain': sc.tl.louvain,
        'leiden': sc.tl.leiden,
    }
    
    if not use_gpu:
        kwargs |= cpu_kwargs
    logging.info(f'{algorithm} clustering with {kwargs}...')
    cluster_func = algorithm_map.get(algorithm, KeyError(f'Unknown clustering algorithm: {algorithm}'))
    cluster_func(adata, **kwargs)
    
    if not max_clusters:
        max_clusters = max(1, int(50 * resolution))
    n_clusters = adata.obs[cluster_key].nunique()
    
    if use_gpu and n_clusters > max_clusters:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        logging.info(
            f'Cluster {cluster_key} has {n_clusters} custers, which is more than {max_clusters}.'
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
overwrite = snakemake.params.get('overwrite', False)
USE_GPU = False

# set parameters for clustering
cluster_key = f'{algorithm}_{resolution}_{level}'
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
    try:
        import subprocess
        assert subprocess.run('nvidia-smi', shell=True).returncode == 0
        from rapids_singlecell.tl import leiden, louvain
        from rapids_singlecell.pp import neighbors
        USE_GPU = True
    except Exception as e:
        logging.info(f'Error importing rapids found...\n{e}')
        from scanpy.tools import leiden, louvain
        from scanpy.preprocessing import neighbors
    
    read_kwargs = dict(obs='obs')
    
    if level <= 1:
        logging.info(f'Read anndata file {input_file}...')
        read_kwargs |= dict(uns='uns', obsp='obsp')
        if input_file.endswith('.h5ad'):
            read_kwargs |= dict(obsm='obsm')
        adata = read_anndata(input_file, **read_kwargs)
        
        # select neighbors
        neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
        logging.info(f'Select neighbors for "{neighbors_key}"...')
        check_and_set_neighbors_key(adata, neighbors_key)
        
        # following observations from https://github.com/rapidsai/cugraph/issues/4072#issuecomment-2074822898
        adata.obsp['connectivities'] = adata.obsp['connectivities'].astype('float64')
        
        adata = apply_clustering(
            adata,
            use_gpu=USE_GPU,
            cpu_kwargs=cpu_kwargs,
            resolution=resolution,
            key_added=cluster_key,
        )
    else:
        logging.info(f'Read anndata file {input_file}...')
        read_kwargs |= dict(obsm='obsm')
        adata = read_anndata(input_file, **read_kwargs)
        
        neighbors_args = dict(use_rep='X_pca') | snakemake.params.get('neighbors_args', {})
        
        # TODO: parallelize using snakemake scripts
        prev_cluster_key = f'{algorithm}_{resolution}_{level-1}'
        for cluster in adata.obs[prev_cluster_key].unique():
            logging.info(f'Subsetting to {prev_cluster_key}={cluster}...')
            sub_adata = adata[adata.obs[prev_cluster_key] == cluster].copy()
            logging.info(sub_adata.shape)
            
            logging.info(f'Compute neighbors for "{neighbors_args}"...')
            neighbors(sub_adata, **neighbors_args) # TODO: custom parameters
            
            logging.info(f'Apply clustering for "{cluster_key}"...')
            sub_adata = apply_clustering(
                sub_adata,
                use_gpu=USE_GPU,
                cpu_kwargs=cpu_kwargs,
                resolution=resolution,
                key_added=cluster_key,
            )
            adata.obs.loc[sub_adata.obs.index, cluster_key] = sub_adata.obs[[prev_cluster_key, cluster_key]].agg('_'.join, axis=1)


logging.info(f'Write {cluster_key} to {output_file}...')
adata.obs = adata.obs[[cluster_key]]
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)