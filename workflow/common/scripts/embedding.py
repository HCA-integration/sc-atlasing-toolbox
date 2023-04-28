import sys
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params

if params is None:
    params = {}

adata = sc.read(input_file)
print(adata)

# plot embedding
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.embedding(adata, **params)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
