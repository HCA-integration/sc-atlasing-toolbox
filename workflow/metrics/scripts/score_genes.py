from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import anndata as ad

from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked
from utils.processing import sc


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
unintegrated_layer = params.get('unintegrated_layer', 'X')
raw_counts_layer = params.get('raw_counts_layer', unintegrated_layer)
gene_set = params['gene_set']
var_key = 'metrics_features'


logger.info(f'Read {input_file} ...')
adata = read_anndata(
    input_file,
    X=unintegrated_layer,
    obs='obs',
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
    adata.obs['score'] = np.nan
    
else:
    # also include genes of interest to subset
    adata.var.loc[adata.var_names.isin(gene_set), var_key] = True

    # subset to HVGs (used for control genes) + genes of interest
    adata, _ = subset_hvg(
        adata,
        var_column=var_key,
        min_cells=1,
        compute_dask=True
    )

    # calculate gene score
    sc.tl.score_genes(adata, gene_list=gene_set)

    # TODO: bin expression (batch-aware?)
    # save to .X

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)