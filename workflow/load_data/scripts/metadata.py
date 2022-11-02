from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import anndata
import scanpy as sc

in_file = snakemake.input.h5ad
out_file = snakemake.output.h5ad
out_plot = snakemake.output.plot
wildcards = snakemake.wildcards
meta = snakemake.params.meta

adata = sc.read(in_file, as_sparse=['X'])
print(adata)

adata.uns['dataset'] = wildcards.dataset
adata.obs['dataset'] = wildcards.dataset
adata.obs['organ'] = meta['organ']
adata.uns['meta'] = meta

donor_column = meta['donor_column']
sample_columns = [s.strip() for s in meta['sample_column'].split('+')]
adata.obs['donor'] = adata.obs[donor_column]
adata.obs['sample'] = adata.obs[sample_columns].apply(lambda x: '-'.join(x), axis=1)

for key, value in meta.items():
    adata.obs[key] = value

# ensure only raw counts are kept
adata.layers['final'] = adata.X.copy()
if isinstance(adata.raw, anndata._core.raw.Raw):
    adata.X = adata.raw.X.copy()
else:
    adata.X = adata.X.copy()
adata.X = csr_matrix(adata.X)

adata.write(out_file, compression='gzip')

# plot count distribution -> save to file
plt.hist(adata.X.data, bins=60)
plt.savefig(out_plot)
