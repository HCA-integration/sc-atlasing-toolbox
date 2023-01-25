import scanpy as sc
from metrics.utils import compute_neighbors, get_from_adata


input_adata = snakemake.input.h5ad
output_file = snakemake.output.h5ad


adata = sc.read(input_adata)
meta = get_from_adata(adata)

for output_type in meta['output_types']:
    # compute neighbors
    ad_knn = compute_neighbors(adata, output_type)
    adata.obsp['connectivities_' + output_type] = ad_knn.obsp['connectivities']
    adata.obsp['distances_' + output_type] = ad_knn.obsp['distances']

print(adata)
adata.write(output_file)