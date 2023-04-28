import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params

if params is None:
    params = {}
else:
    params = {k: v for k,v in params.items()}

print('params:', params)

assert 'basis' in params, '"basis" is positional and must be defined'
basis = params['basis']

adata = sc.read(input_file)
print(adata)


# convert to dense matrix
try:
    adata.obsm[basis] = np.array(adata.obsm[basis].todense())
except:
    pass

# plot embedding
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.embedding(adata, ncols=1, **params)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
