from pathlib import Path
import numpy as np
import pandas as pd
import doubletdetection
import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batch_key = snakemake.params.get('batch_key')
batch = str(snakemake.wildcards.batch)
threads = snakemake.threads

logging.info(f'Read {input_zarr}...')
adata = read_anndata(
    input_zarr,
    X='X',
    obs='obs',
    backed=True,
    dask=True,
)

logging.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key].astype(str) == batch, :].copy()
logging.info(adata.__str__())

if adata.n_obs < 10:
    columns = ['doubletdetection_score', 'doubletdetection_prediction']
    df = pd.DataFrame(index=adata.obs.index, columns=columns, dtype=float).fillna(0)
    df.to_csv(output_tsv, sep='\t')
    exit(0)

# load data to memory
adata = dask_compute(adata)

# run doubletdetection
logging.info('Run doubletdetection...')
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    n_top_var_genes=4000,
    n_components=np.min([adata.n_obs, adata.n_vars, 30]),
    clustering_algorithm="leiden",
    n_jobs=threads,
)
labels = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
scores = clf.doublet_score()

# save results
logging.info('Save results...')
adata.obs['doubletdetection_score'] = scores
adata.obs['doubletdetection_prediction'] = labels
df = adata.obs[['doubletdetection_score', 'doubletdetection_prediction']]
df.to_csv(output_tsv, sep='\t')
