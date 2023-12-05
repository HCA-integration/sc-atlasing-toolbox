from pathlib import Path
import numpy as np
import pandas as pd
import doubletdetection
import logging
logger = logging.getLogger('doubletdetection')

from utils.io import read_anndata, link_zarr

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batch_key = snakemake.params.get('batch_key')
batch = snakemake.wildcards.batch
threads = snakemake.threads

# read data
logger.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, X='X', obs='obs')

# subset to batch
logger.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key] == batch].copy()

if adata.n_obs < 10:
    columns = ['doubletdetection_score', 'doubletdetection_prediction']
    df = pd.DataFrame(index=adata.obs.index, columns=columns, dtype=float).fillna(0)
    df.to_csv(output_tsv, sep='\t')
    exit(0)

# run doubletdetection
logger.info('Run doubletdetection...')
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
logger.info('Save results...')
adata.obs['doubletdetection_score'] = scores
adata.obs['doubletdetection_prediction'] = labels
df = adata.obs[['doubletdetection_score', 'doubletdetection_prediction']]
df.to_csv(output_tsv, sep='\t')
