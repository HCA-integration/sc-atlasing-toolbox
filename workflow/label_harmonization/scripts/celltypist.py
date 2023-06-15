import logging
import scanpy as sc
from scipy import sparse
import celltypist

from utils.io import read_anndata

logger = logging.getLogger(__name__)

input_file = snakemake.input[0]
output_h5ad = snakemake.output.h5ad
output_reannotation = snakemake.output.reannotation
output_relation = snakemake.output.relation
output_model = snakemake.output.model
author_label_key = snakemake.params.author_label_key
dataset_key = snakemake.params.dataset_key
params = snakemake.params.params
subsample = snakemake.params.subsample
force_scale = snakemake.params.force_scale

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
    scaled = adata.uns['preprocessing']['scaled']
except KeyError:
    scaled = False

# scale for PCT if not already scaled
if 'use_pct' in params and params['use_pct'] and (force_scale or not scaled):
    logger.info("Scale .X for PCT")
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

# remove count matrices
del adata.X
del adata.raw
del adata.layers

print(adata.obs)
adata.obs[['low_hierarchy', 'high_hierarchy']] = \
    alignment.reannotation.loc[adata.obs_names, ['reannotation', 'group']]
print(adata.obs)

# remove count matrices
del adata.X
del adata.layers

logger.info('Write file...')
adata.write(output_h5ad, compression='lzf')
