import sys
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = {k: v for k, v in snakemake.params.items()}

adata = sc.read(input_file)
print(adata)

use_rep = 'X'
if 'use_rep' in params.keys():
    use_rep = params['use_rep']
    if use_rep in adata.obsm.keys():
        use_rep = use_rep
    else:
        print(
            f'WARNING: use_rep="{use_rep}" not in adata.obsm, defaulting to use_rep="X"',
            file=sys.stderr
        )
        use_rep = 'X'
    del params['use_rep']

# convert to dense matrix
try:
    if use_rep == 'X':
        adata.X = np.array(adata.X.todense())
    else:
        adata.obsm[use_rep] = np.array(adata.obsm[use_rep].todense())
except:
    print(f'{use_rep} already dense')


# compute neighbors if required
if 'neighbors' not in adata.uns.keys() or adata.uns['neighbors']['params'].get('use_rep') != use_rep:
    print(f'recompute neighbors using {use_rep}...')
    try:
        sc.pp.neighbors(adata, use_rep=use_rep, method='rapids')
    except:
        print('sc.pp.neighbors: Rapids failed, defaulting to UMAP implementation')
        sc.pp.neighbors(adata, use_rep=use_rep)

# compute UMAP
try:
    sc.tl.umap(adata, method='rapids')
except:
    print('sc.tl.umap: Rapids failed, defaulting to UMAP implementation')
    sc.tl.umap(adata)

# plot UMAP
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.umap(adata, **params)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
