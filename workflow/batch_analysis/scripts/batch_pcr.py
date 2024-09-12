import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
import yaml
from concurrent.futures import ProcessPoolExecutor
# from multiprocessing.pool import ThreadPool
try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    print('no hardware acceleration for sklearn', flush=True)

from utils.io import read_anndata


input_file = snakemake.input.anndata
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = snakemake.params.get('sample_key')
n_threads = np.max([snakemake.threads, 1])

print('Read anndata file...', flush=True)
adata = read_anndata(input_file, obsm='obsm', obs='obs', uns='uns')
adata = adata[adata.obs[covariate].notna()].copy()
n_covariate = adata.obs[covariate].nunique()

# set default sample key
if sample_key is None or sample_key == 'None':
    print('Using index as sample key...', flush=True)
    sample_key = 'index'
    adata.obs[sample_key] = adata.obs.index

# make sure the PCA embedding is an array
if not isinstance(adata.obsm['X_pca'], np.ndarray):
    adata.obsm['X_pca'] = adata.obsm['X_pca'].toarray()


print('Read covariate setup...', flush=True)
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
# n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
n_permute = snakemake.params.get('n_permute', 0)
print(f'n_permute: {n_permute}', flush=True)

# PCR for permuted covariates
print(f'Permutate covariate: "{covariate}"...', flush=True)
perm_covariates = []

if covariate == sample_key:
    print('Sample key is the same as covariate, skipping permutation...', flush=True)
    n_permute = 0

for i in range(n_permute):
    # aggregate and permutate covariate
    cov_per_sample = adata.obs.groupby(sample_key, observed=True).agg({covariate: 'first'})
    cov_map = dict(zip(cov_per_sample.index, cov_per_sample[covariate].sample(frac=1)))
    
    # apply permutation
    covariate_perm = f'{covariate}-{i}'
    adata.obs[covariate_perm] = adata.obs[sample_key].map(cov_map)
    perm_covariates.append(covariate_perm)

# defragment adata.obs
adata.obs = adata.obs.copy()

# PC regression for all covariates
covariates = [covariate]+perm_covariates


def compute_pcr(adata, covariate, is_permuted, n_threads=1):
    print(f'PCR for covariate: "{covariate}"', flush=True)
    pcr = scib.me.pcr(
        adata,
        covariate=covariate,
        recompute_pca=False,
        verbose=False,
        linreg_method='numpy',
        n_threads=n_threads,
    )
    print(f'covariate: {covariate}, pcr: {pcr}', flush=True)
    return (covariate, is_permuted, pcr)


with ProcessPoolExecutor(max_workers=n_threads) as executor:
# with ThreadPool(processes=n_threads) as pool:
    chunk_size = max(1, 2 * int(n_permute / n_threads))
    print(f'chunk_size: {chunk_size}', flush=True)
    n_covariates = len(covariates)
    pcr_scores = list(
        executor.map(
            compute_pcr,
            [adata] * n_covariates,
            covariates,
            [x in perm_covariates for x in covariates],
            # [chunk_size] * n_covariates,
            chunksize=chunk_size,
        )
    )

if covariate == sample_key:
    # when covariate is the same as the group variable, save computations and impute same value as covariate score
    perm_score = (
        f'{covariate}-0',
        True,
        pcr_scores[0][2]
    )
    pcr_scores.append(perm_score)

df = pd.DataFrame.from_records(
    pcr_scores,
    columns=['covariate', 'permuted', 'pcr'],
)

# calculate summary stats
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')
df['n_covariates'] = n_covariate
df['perm_mean'] = df.loc[df['permuted'], 'pcr'].mean()
df['perm_std'] = df.loc[df['permuted'], 'pcr'].std()
df['z_score'] = (df['pcr'] - df['perm_mean']) / df['perm_std']
df['signif'] = df['z_score'] > 1.5
df['p-val'] = df.loc[df['permuted'], 'signif'].sum() / n_permute

print(df, flush=True)
df.to_csv(output_file, sep='\t', index=False)