import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')

from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)
neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)

if neighbors_key not in adata.uns:
    assert 'connectivities' in adata.obsp
    assert 'distances' in adata.obsp
    adata.uns[neighbors_key] = {
        'connectivities_key': 'connectivities',
        'distances_key': 'distances',
    }

neighbors = adata.uns[neighbors_key]
adata.uns['neighbors'] = neighbors
adata.obsp['connectivities'] = adata.obsp[neighbors['connectivities_key']]
adata.obsp['distances'] = adata.obsp[neighbors['distances_key']]

logging.info(f'Leiden clustering with resolution {resolution}...')
cluster_key = f'leiden_{resolution}'
sc.tl.leiden(
    adata,
    resolution=resolution,
    key_added=cluster_key,
)

logging.info('Write file...')
adata.obs[cluster_key].to_csv(output_file, sep='\t', index=True)