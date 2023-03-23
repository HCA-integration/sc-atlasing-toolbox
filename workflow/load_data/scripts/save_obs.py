import anndata

in_file = snakemake.input[0]
out_file = snakemake.output[0]

adata = anndata.read_zarr(in_file)
adata.obs.to_csv(out_file, sep='\t', index=False)
