from pathlib import Path
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

if adata.n_obs == 0:
    adata.obs.to_csv(output_tsv, sep='\t')
    exit(0)


# subset to batch
logger.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key] == batch].copy()

# run doubletdetection
logger.info('Run doubletdetection...')
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    n_top_var_genes=4000,
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
