from pathlib import Path
from pprint import pformat
import pandas as pd
import numpy as np
import anndata
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)

from utils_pipeline.io import link_zarr

in_file = snakemake.input[0]
in_dcp = snakemake.input[1]
out_obs = snakemake.output.obs
out_stats = snakemake.output.stats
out_adata = snakemake.output.zarr
metadata_columns = snakemake.params.metadata_cols

def explode_table(df, col, sep=' \|\| '):
    df = df.copy()
    df[col] = df[col].astype(str).str.split(pat=sep)
    return df.explode(col)


id_cols = snakemake.params.id_cols

dcp_tsv = pd.read_table(in_dcp)
with zarr.open(in_file) as z:
    obs_df = read_elem(z['obs'])
n_obs = obs_df.shape[0]

# identify ID column
cols = []
intersection = []
n_donors = []
cxg_ids = []

for col in id_cols:
    if col not in dcp_tsv.columns:
        continue
    
    for cxg_id in ['donor_id', 'sample']:
        # explode ID column
        dcp_tsv_exploded = explode_table(dcp_tsv, col)

        # count ID overlaps
        dcp_ids = set(dcp_tsv_exploded[col].unique())
        obs_ids = set(obs_df[cxg_id].unique())
        intersect = obs_ids.intersection(dcp_ids)

        cols.append(col)
        cxg_ids.append(cxg_id)
        intersection.append(len(intersect))
        n_donors.append(len(obs_ids))

intersect_df = pd.DataFrame(
    {
        'dcp_column': cols,
        'cxg_column': cxg_ids,
        'intersection': intersection,
        'n_cxg': n_donors,
    }
)
intersect_df['intersection_fraction'] = intersect_df['intersection'] / intersect_df['n_cxg']
logging.info(intersect_df)

# identify ID column with most overlap
argmax = intersect_df['intersection'].argmax()
intersect_max = intersect_df.iloc[argmax,:].copy()
assert isinstance(intersect_max, pd.Series), f'max intersection is not a Series\n{intersect_max}'
id_col = intersect_max['dcp_column']
cxg_col = intersect_max['cxg_column']

# get IDs that don't match
intersect_max['mismatched_dcp'] = list(set(dcp_tsv[id_col].unique()) - set(obs_df[cxg_col].unique()))
intersect_max['mismatched_cxg'] = list(set(obs_df[cxg_col].unique()) - set(dcp_tsv[id_col].unique()))
intersect_max['n_mismatched_dcp'] = len(intersect_max['mismatched_dcp'])
intersect_max['n_mismatched_cxg'] = len(intersect_max['mismatched_cxg'])

metadata_columns  = set([id_col] + metadata_columns)
metadata_columns = [c for c in dcp_tsv.columns if sum([c.startswith(mc) for mc in metadata_columns]) > 0]
logging.info('DCP metadata columns:')
logging.info(pformat(metadata_columns))

# explode columns
dcp_tsv = explode_table(dcp_tsv, id_col).dropna(subset=[id_col])
dcp_tsv = dcp_tsv[metadata_columns].drop_duplicates()
logging.info(dcp_tsv)

# deal with columns that cause error when writing zarr file
for col in metadata_columns:
    if dcp_tsv[col].apply(isinstance, args=((str, bool, np.bool_),)).any():
        dcp_tsv[col] = dcp_tsv[col].astype(str)
    elif dcp_tsv[col].isnull().all():
        # remove columns that are all null
        del dcp_tsv[col]

# rename columns
obs_df['index'] = obs_df.reset_index().index

# merge on ID column
obs_df = obs_df.merge(
    dcp_tsv,
    left_on=cxg_col,
    right_on=id_col,
    how='left'
)

# make unique per barcode
obs_df = obs_df.drop_duplicates(subset=["index"])
del obs_df['index']
assert n_obs == obs_df.shape[0], f'Number of observations changed from {n_obs} to {obs_df.shape[0]}'

# save obs
logging.info(obs_df.shape)
logging.info(f'save obs to {out_obs}...')
obs_df.to_csv(out_obs, sep='\t', index=False)

# save stats
intersect_max.to_csv(out_stats, sep='\t', index=True)
# TODO: metadata column completeness
# TODO: aggregatedness of donor ID

# save anndata
logging.info(f'Write to {out_adata}...')
adata = anndata.AnnData(obs=obs_df)
adata.write_zarr(out_adata)

if in_file.endswith('.zarr'):
    input_files = [f.name for f in Path(in_file).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obs']]
    link_zarr(
        in_dir=in_file,
        out_dir=out_adata,
        file_names=files_to_keep,
        overwrite=True,
)
