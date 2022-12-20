"""
Normalisation
"""
import anndata
import scanpy as sc

input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr

adata = anndata.read_zarr(input_zarr)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.write_zarr(output_zarr)
