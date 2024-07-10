from utils.io import read_anndata, write_zarr_linked
from utils.misc import merge


input_clusters = snakemake.input.cluster_anno
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

dfs = [read_anndata(file, obs='obs').obs for file in input_clusters]
cluster_df = merge(dfs, left_index=True, right_index=True)
print(cluster_df, flush=True)

kwargs = dict() if input_zarr.endswith('.h5ad') else dict(obs='obs', uns='uns')
adata = read_anndata(input_zarr, **kwargs)
adata.obs = adata.obs.merge(cluster_df, left_index=True, right_index=True, how='left')
adata.uns['clustering'] = {
    'neighbors_key': snakemake.params.get('neighbors_key', 'neighbors'),
}

write_zarr_linked(
    adata,
    in_dir=input_zarr,
    out_dir=output_zarr,
    files_to_keep=['obs', 'uns'],
)
