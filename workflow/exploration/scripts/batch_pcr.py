import logging
logger = logging.getLogger('Batch PCR')
logger.setLevel(logging.INFO)
import pandas as pd
from scipy import sparse
from matplotlib import pyplot as plt
import seaborn as sns
import scib
import numpy as np

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    logger.info(f'no hardware acceleration for sklearn')

from utils.io import read_anndata

def covariate_invalid(adata, covariate):
    return (covariate not in adata.obs.columns) or (adata.obs[covariate][adata.obs[covariate].notna()].nunique() < 2)

input_file = snakemake.input.zarr
input_metadata = snakemake.input.metadata if 'metadata' in snakemake.input else None
output_barplot = snakemake.output.barplot
dataset = snakemake.params.dataset
covariates = snakemake.params['covariates']
perm_covariates = snakemake.params['permutation_covariates']
n_permute = snakemake.params['n_permute']
sample_key = snakemake.params['sample_key']

adata = read_anndata(input_file)
if input_metadata is not None:
    adata.obs = pd.read_table(input_metadata, low_memory=False)
logger.info(adata.obs.head())

# make sure the PCA embedding is an array
if not isinstance(adata.obsm["X_pca"], np.ndarray):
    adata.obsm["X_pca"] = adata.obsm["X_pca"].toarray()

if adata.n_obs == 0:
    plt.savefig(output_barplot)
    exit(0)

# PCR for permuted covariates
for covariate in perm_covariates:
    logger.info(f'Permute covariate: "{covariate}"')
    if covariate_invalid(adata, covariate):
        logger.info('    covariate not in adata or only 1 unique element, skipping...')
        if covariate in covariates:
            covariates.remove(covariate)
        continue
    # permute n times
    for i in range(n_permute):
        cov_per_sample = adata.obs.groupby(sample_key).agg({covariate: 'first'})
        cov_map = dict(
            zip(
                cov_per_sample.index,
                cov_per_sample[covariate].sample(cov_per_sample.shape[0])
            )
        )
        covariate_perm = f'{covariate}_permuted-{i}'
        adata.obs[covariate_perm] = adata.obs[sample_key].map(cov_map)
        covariates.append(covariate_perm)

# PC regression for all covariates
pcr_scores = []
for covariate in list(covariates):
    logger.info(f'PCR for covariate: "{covariate}"')
    if covariate_invalid(adata, covariate):
        logger.info('    covariate not in adata or only 1 unique element, skipping...')
        if covariate in covariates:
            covariates.remove(covariate)
        continue

    if isinstance(adata.X, (sparse.csr_matrix, sparse.csc_matrix)):
        adata.X = adata.X.todense()

    # pcr = scib_metrics.utils.principal_component_regression(
    #     adata.X,
    #     covariate=adata.obs[covariate],
    #     categorical=True,
    #     n_components=50
    # )

    pcr = scib.me.pcr(
        adata[adata.obs[covariate].notna()],
        covariate=covariate,
        recompute_pca=False,
        verbose=False
    )
    pcr_scores.append(pcr)
    logger.info(pcr)

df = pd.DataFrame.from_dict(dict(covariate=covariates, pcr=pcr_scores))
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')

sns.barplot(
    data=df,
    x='pcr',
    y='covariate',
    dodge=False,
    order=df.sort_values('pcr', ascending=False)['covariate'].unique(),
).set(title=f'PCR of covariates for: {dataset}')
plt.xticks(rotation=90)

plt.savefig(output_barplot, bbox_inches='tight')
