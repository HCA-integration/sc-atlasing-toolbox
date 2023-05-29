from matplotlib import pyplot as plt
import scanpy as sc

input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params

if params is None:
    params = {}
else:
    params = {k: v for k,v in params.items()}

adata = sc.read(input_file)

sc.pl.dotplot(
    adata,
    **params,
    show=False
)
plt.savefig(output_file, bbox_inches='tight', dpi=200)