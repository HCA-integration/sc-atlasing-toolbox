import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
from anndata.experimental import read_elem
import zarr
import yaml
from multiprocessing.pool import ThreadPool

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    logger.info('no hardware acceleration for sklearn')

import logging
logger = logging.getLogger('Batch PCR')
logger.setLevel(logging.INFO)

input_file = snakemake.input.anndata
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = snakemake.params.get('sample_key')
n_threads = snakemake.threads

logger.info('Read anndata file...')
z = zarr.open(input_file)
X_pca = read_elem(z['obsm/X_pca'])
obs = read_elem(z['obs'])
uns = read_elem(z['uns'])
n_covariate = obs[covariate].nunique()

# set default sample key
if sample_key is None or sample_key == 'None':
    sample_key = 'index'
    obs[sample_key] = obs.index

# make sure the PCA embedding is an array
if not isinstance(X_pca, np.ndarray):
    X_pca = X_pca.toarray()


logger.info('Read covariate setup...')
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
logger.info(f'n_permute: {n_permute}')


# PCR for permuted covariates
logger.info(f'Permute covariate: "{covariate}"')
perm_covariates = []
for i in range(n_permute):
    covariate_perm = f'{covariate}-{i}'
    cov_per_sample = obs.groupby(sample_key).agg({covariate: 'first'})
    cov_map = dict(
        zip(
            cov_per_sample.index,
            cov_per_sample[covariate].sample(cov_per_sample.shape[0])
        )
    )
    # permute covariate
    obs[covariate_perm] = obs[sample_key].map(cov_map)
    perm_covariates.append(covariate_perm)


# PC regression for all covariates
adata = AnnData(obs=obs, obsm={'X_pca': X_pca}, uns=uns)
covariates = [covariate]+perm_covariates
# adatas = [adata[obs[covariate].notna()] for _ in covariates]
pcr_scores = []


def compute_pcr(covariate):
    logger.info(f'PCR for covariate: "{covariate}"')
    pcr = scib.me.pcr(
        adata[obs[covariate].notna()],
        covariate=covariate,
        recompute_pca=False,
        verbose=False
    )
    logger.info(pcr)
    return pcr


with ThreadPool(processes=n_threads) as pool:
    for pcr in pool.map(compute_pcr, covariates):
        pcr_scores.append(pcr)

df = pd.DataFrame.from_dict(
    {
        'covariate': covariates,
        'pcr': pcr_scores,
        'permuted': [c in perm_covariates for c in covariates],
    }
)
df['n_covariates'] = f'n={n_covariate}'
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')

df.to_csv(output_file, sep='\t', index=False)