import anndata


in_file = snakemake.input.zarr
out_file = snakemake.output.h5ad

anndata.read_zarr(in_file).write_h5ad(out_file)
