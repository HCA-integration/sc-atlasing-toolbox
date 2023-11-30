from pathlib import Path
import pandas as pd
import scanpy as sc
import logging
logger = logging.getLogger('Scrublet')

from utils.io import read_anndata, link_zarr

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batch_key = snakemake.params.get('batch_key')
batch = snakemake.wildcards.batch

# read data
logger.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, X='X', obs='obs')

if adata.n_obs == 0:
    adata.obs.to_csv(output_tsv, sep='\t')
    exit(0)

# subset to batch
logger.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key] == batch].copy()

# run scrublet
logger.info('Run scrublet...')
sc.external.pp.scrublet(
    adata,
    batch_key=None,
    sim_doublet_ratio=2.0,
    expected_doublet_rate=0.05,
    stdev_doublet_rate=0.02,
    synthetic_doublet_umi_subsampling=1.0,
    knn_dist_metric='euclidean',
    normalize_variance=True,
    log_transform=False,
    mean_center=True,
    n_prin_comps=30,
    use_approx_neighbors=True,
    get_doublet_neighbor_parents=False,
    n_neighbors=None,
    threshold=None,
    verbose=True,
    copy=False,
    random_state=0,
)

# save results
logger.info('Save results...')
df = adata.obs[['doublet_score', 'predicted_doublet']].rename(
    columns={
        'doublet_score': 'scrublet_score',
        'predicted_doublet': 'scrublet_prediction',
    }
)
print(df)
df.to_csv(output_tsv, sep='\t')
