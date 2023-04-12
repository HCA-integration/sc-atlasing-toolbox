from scipy import sparse
import celltypist

from utils.io import read_anndata


input_file = snakemake.input[0]
author_label_key = snakemake.params.author_label_key
dataset_key = snakemake.params.dataset_key
params = snakemake.params.params
output_h5ad = snakemake.output.h5ad
output_reannotation = snakemake.output.reannotation
output_relation = snakemake.output.relation
output_model = snakemake.output.model

print('read...')
adata = read_anndata(input_file)

# subset to HVGs
if 'highly_variable' in adata.var:
    adata = adata[:, adata.var['highly_variable']].copy()

# if sparse.issparse(adata.X):
#     # adata.X.indices = adata.X.indices.astype('int32')
#     # adata.X.indptr = adata.X.indptr.astype('int32')
#     adata.X = adata.X.toarray()

print(adata)

alignment = celltypist.harmonize(
    adata,
    dataset_key,
    author_label_key,
    reannotate=True,
    **params
)
alignment.write(output_model)
alignment.reannotation.to_csv(output_reannotation, sep='\t')

df = alignment.relation
df['group'] = alignment.groups
print(df)
df.to_csv(output_relation, sep='\t', index=False)

print(adata.obs)
adata.obs = adata.obs.merge(
    alignment.reannotation[['reannotation', 'group']],
    left_index=True,
    right_index=True,
    how='left'
)
print(adata.obs)
adata.write(output_h5ad)