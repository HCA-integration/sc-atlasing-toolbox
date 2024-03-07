import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
import yaml
from multiprocessing.pool import ThreadPool
import logging
logger = logging.getLogger('Batch PCR')

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    logger.error('no hardware acceleration for sklearn')

from utils.io import read_anndata


input_file = snakemake.input.anndata
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = snakemake.params.get('sample')
n_threads = np.max([snakemake.threads, 1])

logger.error('Read anndata file...')
adata = read_anndata(input_file, obsm='obsm', obs='obs', uns='uns')
n_covariate = adata.obs[covariate].nunique()

# set default sample key
if sample_key is None or sample_key == 'None':
    logger.info('Using index as sample key...')
    sample_key = 'index'
    adata.obs[sample_key] = adata.obs.index

# make sure the PCA embedding is an array
if not isinstance(adata.obsm['X_pca'], np.ndarray):
    adata.obsm['X_pca'] = adata.obsm['X_pca'].toarray()


logger.error('Read covariate setup...')
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
# n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
n_permute = snakemake.params.get('n_permute', 0)
logger.error(f'n_permute: {n_permute}')


# PCR for permuted covariates
logger.error(f'Permute covariate: "{covariate}"')
perm_covariates = []
for i in range(n_permute):
    covariate_perm = f'{covariate}-{i}'
    cov_per_sample = adata.obs.groupby(sample_key, observed=True).agg({covariate: 'first'})
    cov_map = dict(
        zip(
            cov_per_sample.index,
            cov_per_sample[covariate].sample(cov_per_sample.shape[0])
        )
    )
    # permute covariate
    adata.obs[covariate_perm] = adata.obs[sample_key].map(cov_map)
    perm_covariates.append(covariate_perm)


# PC regression for all covariates
# adata = AnnData(obs=obs, obsm={'X_pca': X_pca}, uns=uns)
covariates = [covariate]+perm_covariates
# adatas = [adata[obs[covariate].notna()] for covariate in covariates]
pcr_scores = []


def compute_pcr(covariate):
    logger.error(f'PCR for covariate: "{covariate}"')
    pcr = scib.me.pcr(
        adata[adata.obs[covariate].notna()],
        covariate=covariate,
        recompute_pca=False,
        verbose=False
    )
    logger.error(pcr)
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