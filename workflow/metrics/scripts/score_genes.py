from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import anndata as ad
from tqdm import tqdm
from scipy import sparse

from metrics.bootstrap import get_bootstrap_adata
from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked
from utils.processing import sc
from utils.misc import dask_compute


def bins_by_quantiles(matrix: [np.array, np.matrix], n_quantiles: int):
    assert n_quantiles > 1
    
    def bin_values(x):
        bins = np.linspace(0, 1, n_quantiles + 1)
        quantiles = np.quantile(x, bins)
        binned_idx = np.digitize(x, quantiles) - 1
        return bins[binned_idx]
    
    return np.apply_along_axis(bin_values, axis=0, arr=matrix)


input_file = snakemake.input[0]
output_file = snakemake.output[0]
batch_key = snakemake.wildcards.batch
gene_set_name = snakemake.wildcards.gene_set

params = snakemake.params
unintegrated_layer = params.get('unintegrated_layer', 'X')
raw_counts_layer = params.get('raw_counts_layer', unintegrated_layer)
gene_set = params['gene_set']
n_quantiles = params.get('n_quantiles', 2)
var_key = 'metrics_features'

files_to_keep = ['obs', 'obsm']   # TODO: only overwrite specific obsm slots


logger.info(f'Read {input_file} ...')
adata = read_anndata(
    input_file,
    X=unintegrated_layer,
    obs='obs',
    obsm='obsm',
    var='var',
    dask=True,
    backed=True,
)

adata.layers['raw_counts'] = read_anndata(
    input_file,
    X=raw_counts_layer,
    dask=True,
    backed=True,
).X
 
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']

# filter genes of interest to genes in dataset
gene_set = [g for g in gene_set if g in adata.var_names]

if len(gene_set) == 0:
    logging.warning('No genes of interest found in dataset')
    adata.obs[gene_set_name] = np.nan
    adata.obsm['random_gene_scores'] = np.empty((adata.n_obs, 0))
    adata.obsm['binned_expression'] = sparse.csr_matrix((adata.n_obs, 0))
    
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
        files_to_keep=files_to_keep,
    )
    exit(0)


# subset adata if dataset too large TODO: move to prepare script?
n_subset = int(4e6)
if adata.n_obs > n_subset:
    adata.obsp = read_anndata(input_file, obs='obs', obsp='obsp').obsp
    files_to_keep.extend(['obsp', 'X', 'layers'])
    adata = get_bootstrap_adata(adata, size=n_subset)

# subset to HVGs (used for control genes) + genes of interest
adata.var.loc[adata.var_names.isin(gene_set), var_key] = True
adata = adata[:, adata.var[var_key]].copy()
adata = dask_compute(adata, layers='X')

logging.info(f'Gene score for gene set with {len(gene_set)} genes...')
sc.tl.score_genes(
    adata,
    gene_list=gene_set,
    score_name=gene_set_name
)

logging.info('Add random gene scores...')
scores = []
n_samples = 5
for _ in tqdm(range(n_samples)):
    random_genes = np.random.choice(
        adata.var_names,
        size=len(gene_set),
        replace=False
    )
    sc.tl.score_genes(adata, gene_list=random_genes)
    scores.append(adata.obs['score'].values)
adata.obsm['random_gene_scores'] = np.array(scores).T
del adata.obs['score']

logging.info('Bin expression...')
adata = dask_compute(adata[:, gene_set].copy(), layers='raw_counts')
adata.obsm[f'binned_expression'] = bins_by_quantiles(
    adata.layers['raw_counts'].toarray(),
    n_quantiles=n_quantiles,
)

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
)