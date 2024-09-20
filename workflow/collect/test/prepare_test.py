import anndata as ad

file_path = 'input/pbmc68k.h5ad'
out_path = 'input/pbmc68k_modified.zarr'

print(f'Read {file_path}...', flush=True)
adata = ad.read_h5ad(file_path)
adata.obs['louvain'] = adata.obs['batch']
adata.obs['new_column'] = 'for collect module'

print(f'Write to {out_path}...', flush=True)
adata.write_zarr(out_path)