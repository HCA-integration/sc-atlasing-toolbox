from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

import pandas as pd
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import anndata
import scanpy as sc
import numpy as np
from utils import SCHEMAS, get_union

in_file = snakemake.input.h5ad
out_file = snakemake.output.zarr
out_plot = snakemake.output.plot
wildcards = snakemake.wildcards
meta = snakemake.params.meta
schema_file = snakemake.input.schema

logging.info('meta:')
logging.info(pformat(meta))

# optional annotations file
annotation_file = meta['annotation_file']
if not pd.isna(annotation_file):
    logging.info(f'check if annotation file {annotation_file} exists')
    assert Path(annotation_file).exists()

# h5ad
logging.info(f'\033[0;36mread\033[0m {in_file}...')
try:
    adata = sc.read(in_file, as_sparse=['X'])
except:
    adata = sc.read_loom(in_file, sparse=True)
logging.info(adata)

# Adding general dataset info to uns and obs
adata.uns['meta'] = meta
for meta_i in ["organ", "study", "dataset"]:
    adata.obs[meta_i] = meta[meta_i]
    adata.uns[meta_i] = meta[meta_i]

# add annotation if available
author_annotation = meta['author_annotation']
if not pd.isna(annotation_file):
    logging.info(f'Add annotations from {annotation_file}...')
    annotation = pd.read_csv(annotation_file)
    annotation.index = annotation[meta['barcode_column']].astype(str)
    adata.obs.index = adata.obs.index.astype(str)
    # remove column if existing to avoid conflict
    if author_annotation in adata.obs.columns:
        logging.info(f'column {author_annotation} already exists, removing it')
        del adata.obs[author_annotation]
    adata.obs = adata.obs.join(annotation[author_annotation], how='left')
    logging.info(adata.obs[author_annotation])

# assign sample and donor variables
donor_column = meta['donor_column']
adata.obs['donor'] = adata.obs[donor_column]

sample_columns = [s.strip() for s in meta['sample_column'].split('+')]
adata.obs['sample'] = adata.obs[sample_columns].apply(lambda x: '-'.join(x), axis=1)

# CELLxGENE specific
if 'batch_condition' in adata.uns.keys():
    batch_columns = adata.uns['batch_condition']
    adata.obs['batch'] = adata.obs[batch_columns].apply(lambda x: '-'.join(x), axis=1)
else:
    adata.obs['batch'] = meta['study']

# Checking schema version
if 'schema_version' not in adata.uns.keys():
    adata.uns['schema_version'] = '0.0.0'
if adata.uns['schema_version'] == '2.0.0':
    adata.obs['self_reported_ethnicity'] = adata.obs['ethnicity']
    adata.obs['self_reported_ethnicity_ontology_term_id'] = adata.obs['ethnicity_ontology_term_id']
    adata.obs['donor_id'] = adata.obs['donor']

# Assigning other keys in meta to obs
for key, value in meta.items():
    if isinstance(value, list):
        continue
    adata.obs[key] = value

# ensure raw counts are kept in .X
adata.layers['final'] = adata.X
if isinstance(adata.raw, anndata._core.raw.Raw):
    adata.X = adata.raw.X
else:
    adata.X = adata.X
adata.X = csr_matrix(adata.X)

# add author annotations column
adata.obs['author_annotation'] = adata.obs[author_annotation]
# use author annotations if no cell ontology available
if 'cell_type' not in adata.obs.columns:
    adata.obs['cell_type'] = 'nan'
if adata.obs['cell_type'].nunique() == 1:
    adata.obs['cell_type'] = adata.obs['author_annotation']

# save barcodes in separate column
adata.obs['barcode'] = adata.obs_names
adata.obs_names = adata.uns['dataset'] + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# schemas translation
schemas_df = pd.read_table(schema_file).dropna()
logging.info(schemas_df)
from_schema = meta['schema']
assert from_schema in schemas_df.columns
to_schema = 'cellxgene'
SCHEMAS["NAMES"] = dict(zip(schemas_df[from_schema], schemas_df[to_schema]))
adata.obs.rename(SCHEMAS["NAMES"], inplace=True)

# making sure all columns are in the object
all_columns = get_union(SCHEMAS["CELLxGENE_OBS"], SCHEMAS["EXTRA_COLUMNS"])
for column in all_columns:
    if column not in adata.obs.columns:
        adata.obs[column] = np.nan
# keep only relevant columns
adata.obs = adata.obs[all_columns].copy()

# make sure all vars are present
if "feature_name" not in adata.var.columns:
    adata.var["feature_name"] = adata.var_names.tolist()
for column in SCHEMAS["CELLxGENE_VARS"]:
    if column not in adata.var.columns:
        adata.var[column] = np.nan
adata.var = adata.var[SCHEMAS["CELLxGENE_VARS"]]

if 'feature_id' not in adata.var.columns:
    adata.var['feature_id'] = adata.var_names
adata.var.index.set_names('feature_id', inplace=True)

adata.write_zarr(out_file)

# plot count distribution -> save to file
plt.hist(adata.X.data, bins=60)
plt.savefig(out_plot)
