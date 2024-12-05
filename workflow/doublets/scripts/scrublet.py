from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.processing import sc, USE_GPU
from utils.misc import dask_compute

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
layer = snakemake.params.get('layer', 'X')
batch_key = snakemake.params.get('batch_key')
batch = str(snakemake.wildcards.batch)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    backed=True,
    dask=True,
)

logging.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key].astype(str) == batch, :].copy()
logging.info(adata.__str__())

if adata.n_obs < 10:
    columns = ['scrublet_score', 'scrublet_prediction']
    df = pd.DataFrame(index=adata.obs.index, columns=columns, dtype=float).fillna(0)
    df.to_csv(output_tsv, sep='\t')
    exit(0)

# load data to memory
adata = dask_compute(adata)

# run scrublet
logging.info('Run scrublet...')

if USE_GPU:
    sc.get.anndata_to_GPU(adata)

sc.pp.scrublet(
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
    n_prin_comps=np.min([adata.n_obs-1, adata.n_vars-1, 30]),
    use_approx_neighbors=True,
    get_doublet_neighbor_parents=False,
    n_neighbors=None,
    threshold=None,
    verbose=True,
    copy=False,
    random_state=0,
)

# save results
logging.info('Save results...')
df = adata.obs[['doublet_score', 'predicted_doublet']].rename(
    columns={
        'doublet_score': 'scrublet_score',
        'predicted_doublet': 'scrublet_prediction',
    }
)
print(df)
df.to_csv(output_tsv, sep='\t')
