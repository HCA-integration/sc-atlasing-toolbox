from pathlib import Path
import numpy as np
import pandas as pd
import doubletdetection
import anndata as ad
import logging
logger = logging.getLogger('doubletdetection')

from utils.io import read_anndata

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batch_key = snakemake.params.get('batch_key')
batch = str(snakemake.wildcards.batch)
threads = snakemake.threads

logger.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, backed=True, X='X', obs='obs', var='var')

logger.info(f'Subset to batch {batch}...')
if batch_key in adata.obs.columns:
    adata = adata[adata.obs[batch_key].astype(str) == batch, :]
else:
    adata = adata

if isinstance(adata.X, (ad.experimental.CSRDataset, ad.experimental.CSCDataset)):
    adata.X = adata.X.to_memory()

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
