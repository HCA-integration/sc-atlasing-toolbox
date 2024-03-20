import pandas as pd

from utils.io import read_anndata, link_zarr_partial
from utils.misc import merge


input_tsv = snakemake.input.tsv
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

dfs = [pd.read_table(file, index_col=0, dtype='str') for file in input_tsv]
cluster_df = merge(dfs, left_index=True, right_index=True)
print(cluster_df)

adata = read_anndata(input_zarr, obs='obs', uns='uns')
adata.obs = adata.obs.merge(cluster_df, left_index=True, right_index=True, how='left')
adata.uns['clustering'] = {
    'neighbors_key': snakemake.params.get('neighbors_key', 'neighbors'),
    'algorithm': snakemake.params.get('algorithm', 'louvain'),
}

adata.write_zarr(output_zarr)
link_zarr_partial(input_zarr, output_zarr, files_to_keep=['obs', 'uns'])