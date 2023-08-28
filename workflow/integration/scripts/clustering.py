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
from metrics.utils import get_from_adata, select_neighbors

input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file)
meta = get_from_adata(adata)

cluster_keys = []
for output_type in meta['output_types']:
    logging.info(f'Leiden clustering with resolution {resolution} and output type {output_type}...')
    cluster_key = f'leiden_{output_type}_{resolution}'
    cluster_keys.append(cluster_key)
    adata = select_neighbors(adata, output_type)
    sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=cluster_key,
        )

logging.info('Write file...')
adata.obs[cluster_keys].to_csv(output_file, sep='\t', index=True)

if 'cp' in globals():
    cp.cuda.Stream.null.synchronize()
    logging.info('Synchronised parallel jobs')

# hack: manually terminate job
# jid = snakemake.params.get('SLURM_JOB_ID')
# if jid is not None:
#     run(['scancel', jid], check=False)