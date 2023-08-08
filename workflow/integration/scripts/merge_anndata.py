import anndata as ad
import mudata as mu

from utils.io import read_anndata


input_files = snakemake.input
output_file = snakemake.output.h5mu
lineages = snakemake.input.keys()

print(lineages)

adatas = [read_anndata(file) for file in input_files]
mudata = mu.MuData({str(lineage): adata for lineage, adata in zip(lineages, adatas)})
mudata.uns = adatas[0].uns

mudata.write(output_file)
