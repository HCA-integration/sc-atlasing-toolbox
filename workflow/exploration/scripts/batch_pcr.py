import pandas as pd
from scipy import sparse
from matplotlib import pyplot as plt
import seaborn as sns
import anndata
import scib

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    print(f'no hardware acceleration for sklearn')

from utils.io import read_anndata

input_file = snakemake.input.zarr
input_metadata = snakemake.input.metadata
output_barplot = snakemake.output.barplot
dataset = snakemake.params.dataset
covariates = snakemake.params['covariates']
perm_covariates = snakemake.params['permutation_covariates']
sample_key = snakemake.params['sample_key']

adata = read_anndata(input_file)
obs = pd.read_table(input_metadata)

adata.obs = obs

if adata.n_obs == 0:
    plt.savefig(output_barplot)
    exit(0)

# PCR for permuted covariates
for covariate in perm_covariates:
    if adata.obs[covariate].nunique() == 1:
        continue
    # permute 10 times
    for i in range(10):
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
    if adata.obs[covariate].nunique() == 1:
        print(f'skip {covariate}')
        covariates.remove(covariate)
        continue
    # X = adata.X
    # if isinstance(adata.X, (sparse.csr_matrix, sparse.csc_matrix)):
    #     X = X.todense()
    # pcr = scib_metrics.utils.principal_component_regression(
    #     X,
    #     covariate=adata.obs[covariate],
    #     categorical=True,
    #     n_components=50
    # )
    pcr = scib.me.pcr(adata, covariate=covariate, recompute_pca=False, verbose=False)
    pcr_scores.append(pcr)
    print(covariate, pcr)

df = pd.DataFrame.from_dict(dict(covariate=covariates, pcr=pcr_scores))
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')

sns.barplot(
    data=df,
    x='covariate',
    y='pcr',
    dodge=False,
    order=df.sort_values('pcr', ascending=False)['covariate'].unique(),
).set(title=f'PCR of covariates for: {dataset}')
plt.xticks(rotation=90)

plt.savefig(output_barplot, bbox_inches='tight')
