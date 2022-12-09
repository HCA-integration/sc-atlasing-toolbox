from pathlib import Path

import pandas as pd
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import anndata
import scanpy as sc
from utils import CELLxGENE_OBS, CELLxGENE_VARS, EXTRA_COLUMNS, get_union

in_file = snakemake.input.h5ad
out_file = snakemake.output.h5ad
out_plot = snakemake.output.plot
wildcards = snakemake.wildcards
meta = snakemake.params.meta
print(meta)

# optional annotations file
annotation_file = meta['annotation_file']
if not pd.isna(annotation_file):
    print(f'check if annotation file {annotation_file} exists')
    assert Path(annotation_file).exists()

# h5ad
print(f'read {in_file}...')
adata = sc.read(in_file, as_sparse=['X'])
print(adata)

adata.uns['meta'] = meta
adata.uns['dataset'] = meta['dataset']
adata.uns['organ'] = meta['organ']
adata.uns['dataset'] = meta['dataset']

adata.obs['organ'] = meta['organ']
adata.obs['study'] = meta['study']
adata.obs['dataset'] = meta['dataset']

# add annotation if available
author_annotation = meta['author_annotation']
if not pd.isna(annotation_file):
    print(f'Add annotations from {annotation_file}...')
    annotation = pd.read_csv(annotation_file)
    annotation.index = annotation[meta['barcode_column']]
    # remove column if existing to avoid conflict
    if author_annotation in adata.obs.columns:
        del adata.obs[author_annotation]
    adata.obs = adata.obs.join(annotation[author_annotation], how='left')
    print(adata.obs)

# assign sample and donor variables
donor_column = meta['donor_column']
sample_columns = [s.strip() for s in meta['sample_column'].split('+')]
adata.obs['donor'] = adata.obs[donor_column]
adata.obs['sample'] = adata.obs[sample_columns].apply(lambda x: '-'.join(x), axis=1)

if adata.uns['schema_version'] == '2.0.0':
    adata.obs['self_reported_ethnicity'] = adata.obs['ethnicity']
    adata.obs['self_reported_ethnicity_ontology_term_id'] = adata.obs['ethnicity_ontology_term_id']
    adata.obs['donor_id'] = adata.obs['donor']

for key, value in meta.items():
    adata.obs[key] = value

# ensure only raw counts are kept
adata.layers['final'] = adata.X.copy()
if isinstance(adata.raw, anndata._core.raw.Raw):
    adata.X = adata.raw.X.copy()
else:
    adata.X = adata.X.copy()
adata.X = csr_matrix(adata.X)

# add author annotations
adata.obs['author_annotation'] = adata.obs[author_annotation]
# use author annotations if no cell ontology available
if adata.obs['cell_type'].nunique() == 1:
    adata.obs['cell_type'] = adata.obs['author_annotation']

# save barcodes in separate column
adata.obs['barcode'] = adata.obs_names
adata.obs_names = adata.uns['dataset'] + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# keep only relevant columns
adata.obs = adata.obs[get_union(CELLxGENE_OBS, EXTRA_COLUMNS)].copy()
adata.var = adata.var[CELLxGENE_VARS]
adata.var.index.set_names('feature_id', inplace=True)

adata.write(out_file, compression='gzip')

# plot count distribution -> save to file
plt.hist(adata.X.data, bins=60)
plt.savefig(out_plot)
