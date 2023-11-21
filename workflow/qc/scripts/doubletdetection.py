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

# read data
logger.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, X='X', obs='obs')

# subset to batch
logger.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key] == batch].copy()

# run doubletdetection
logger.info('Run doubletdetection...')
clf = doubletdetection.BoostClassifier()
labels = clf.fit(adata.X).predict()
scores = clf.doublet_score()

# save results
logger.info('Save results...')
adata.obs['doubletdetection_score'] = scores
adata.obs['doubletdetection_prediction'] = labels
df = adata.obs[['doubletdetection_score', 'doubletdetection_prediction']]
df.to_csv(output_tsv, sep='\t')
