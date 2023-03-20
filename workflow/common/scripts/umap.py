import sys
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc

input_file = snakemake.input[0]
output_file = snakemake.output[0]

adata = sc.read(input_file)
print(adata)

use_rep = 'X'
if 'use_rep' in snakemake.params.keys():
    use_rep = snakemake.params['use_rep']
    if use_rep in adata.obsm.keys():
        use_rep = use_rep
    else:
        print(
            f'WARNING: use_rep="{use_rep}" not in adata.obsm, defaulting to use_rep="X"',
            file=sys.stderr
        )
        use_rep = 'X'

# convert to dense matrix
try:
    if use_rep == 'X':
        adata.X = np.array(adata.X.todense())
    else:
        adata.obsm[use_rep] = np.array(adata.obsm[use_rep].todense())
except:
    print(f'{use_rep} already dense')


# compute neighbors if required
if 'neighbors' not in adata.uns.keys():
    print(f'recompute neighbors using {use_rep}...')
    try:
        sc.pp.neighbors(adata, use_rep=use_rep, method='rapids')
    except:
        print('Rapids failed, defaulting to UMAP implementation')
        sc.pp.neighbors(adata, use_rep=use_rep)

# compute UMAP
try:
    sc.tl.umap(adata, method='rapids')
except:
    print('Rapids failed, defaulting to UMAP implementation')
    sc.tl.umap(adata)

# pick colours for plot
if 'color' in snakemake.params.keys():
    color = snakemake.params['color']
else:
    color = None

# plot UMAP
sc.set_figure_params(frameon=False, vector_friendly=True, fontsize=9)
sc.pl.umap(adata, color=color, ncols=1)
plt.savefig(output_file, bbox_inches='tight', dpi=200)
