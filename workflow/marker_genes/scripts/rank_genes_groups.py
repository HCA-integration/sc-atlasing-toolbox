import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
rank_gene_func = sc.tl.rank_genes_groups

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
# from utils.processing import sc, USE_GPU

# if USE_GPU:
#     rank_gene_func = sc.tl.rank_genes_groups_logreg
# else:
#     rank_gene_func = sc.tl.rank_genes_groups


input_file = snakemake.input[0]
output_file = snakemake.output.zarr
output_tsv = snakemake.output.tsv
group_key = snakemake.wildcards.group
args = snakemake.params.get('args', {})

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

adata, _ = subset_hvg(adata, var_column='highly_variable')

# filter groups
groups_counts = adata.obs[group_key].value_counts()
groups = groups_counts[groups_counts > 1].index.tolist()

if 'feature_names' in adata.var.columns:
    adata.var_names = adata.var['feature_names']

logging.info(f'Running marker genes analysis for {group_key} and args={args}...')
key = f'marker_genes_group={group_key}'
rank_gene_func(
    adata,
    groupby=group_key,
    groups=groups[0],
    pts=True,
    use_raw=False,
    key_added=key,
    **args
)

logging.info(f'Writing {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns'],
)


logging.info('Create DEG dataframe...')
result = adata.uns[key]
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