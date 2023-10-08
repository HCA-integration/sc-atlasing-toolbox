import numpy as np
import pickle as pkl
import pandas as pd
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata

# try:
#     import rapids_singlecell as rsc
#     logging.info('Using rapids_singlecell')
#     rank_gene_func = rsc.tl.rank_genes_groups_logreg
# except ImportError:
logging.info('Using scanpy')
import scanpy as sc
rank_gene_func = sc.tl.rank_genes_groups


input_file = snakemake.input[0]
input_clusters = snakemake.input.clusters
output_tsv = snakemake.output.tsv
output_uns = snakemake.output.uns
output_rankplot = snakemake.output.rankplot
resolution = snakemake.wildcards.resolution
output_type = snakemake.wildcards.output_type


logging.info('Reading data...')
adata = read_anndata(input_file)
clusters_df = pd.read_table(input_clusters, sep='\t', index_col=0)
cluster_key = clusters_df.columns[0]
adata.obs[cluster_key] = clusters_df[cluster_key].astype(str).astype('category')

# filter groups
groups_counts = adata.obs[cluster_key].value_counts()
groups = groups_counts[groups_counts > 1].index.tolist()


if 'feature_names' in adata.var.columns:
    adata.var_names = adata.var['feature_names']

logging.info(f'Running marker genes analysis for {resolution}...')
rank_gene_func(
    adata,
    groupby=cluster_key,
    groups=groups,
    method='wilcoxon',
    n_genes=200,
    pts=True,
    use_raw=False,
)
result = adata.uns['rank_genes_groups']

with open(output_uns, 'wb') as f:
    pkl.dump(result, f)

logging.info('Create dataframe...')
markers_df = []
for cluster in result['names'].dtype.names:
    current = pd.DataFrame(
        {
            "gene": result["names"][cluster],
            "z-score": result["scores"][cluster],
            "logfoldchange": result["logfoldchanges"][cluster],
            "pval": result["pvals"][cluster],
            "pval_adj": result["pvals_adj"][cluster], 
            "pct_within": result["pts"].loc[result["names"][cluster]][cluster],
            "pct_outside": result["pts_rest"].loc[result["names"][cluster]][cluster],
            "cluster": cluster
        }
    )
    markers_df.append(current)
markers_df = pd.concat(markers_df).set_index('cluster')
markers_df['-log10 pvalue'] = -np.log10(markers_df['pval'])

markers_df.to_csv(output_tsv, sep='\t')

logging.info('Plot...')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig(output_rankplot, dpi=80, bbox_inches='tight')
