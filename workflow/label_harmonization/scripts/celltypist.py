import scanpy as sc
from scipy import sparse
import celltypist

from utils.io import read_anndata


input_file = snakemake.input[0]
output_h5ad = snakemake.output.h5ad
output_reannotation = snakemake.output.reannotation
output_relation = snakemake.output.relation
output_model = snakemake.output.model
author_label_key = snakemake.params.author_label_key
dataset_key = snakemake.params.dataset_key
params = snakemake.params.params
subsample = snakemake.params.subsample

print(params)

print('read...')
adata = read_anndata(input_file)

# subsample genes
if subsample:
    assert 0 < subsample < 1
    sc.pp.subsample(adata, fraction=subsample)

# subset to HVGs
if 'highly_variable' in adata.var:
    adata = adata[:, adata.var['highly_variable']].copy()

try:
    scale = adata.uns['preprocessing']['scaled']
except KeyError:
    scale = False

# scale for PCT if not already scaled
if 'use_pct' in params and not scale:
    sc.pp.scale(adata, max_value=10)

print(adata)

alignment = celltypist.harmonize(
    adata=adata,
    dataset=dataset_key,
    cell_type=author_label_key,
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