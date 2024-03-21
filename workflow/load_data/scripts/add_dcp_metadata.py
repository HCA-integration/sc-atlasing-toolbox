from pathlib import Path
from pprint import pformat
import pandas as pd
import numpy as np
import anndata
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
in_dcp = snakemake.input.dcp
output_file = snakemake.output.zarr
out_stats = snakemake.output.stats
metadata_columns = snakemake.params.metadata_cols

def explode_table(df, col, sep=' \|\| '):
    df = df.copy()
    df[col] = df[col].astype(str).str.split(pat=sep)
    return df.explode(col)


id_cols = snakemake.params.id_cols

dcp_tsv = pd.read_table(in_dcp)
adata = read_anndata(input_file, obs='obs')

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
        obs_ids = set(adata.obs[cxg_id].unique())
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
intersect_max['mismatched_dcp'] = list(set(dcp_tsv[id_col].unique()) - set(adata.obs[cxg_col].unique()))
intersect_max['mismatched_cxg'] = list(set(adata.obs[cxg_col].unique()) - set(dcp_tsv[id_col].unique()))
intersect_max['n_mismatched_dcp'] = len(intersect_max['mismatched_dcp'])
intersect_max['n_mismatched_cxg'] = len(intersect_max['mismatched_cxg'])

if intersect_max['intersection'] > 0:
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
    adata.obs['index'] = adata.obs.reset_index().index

    # merge on ID column
    adata.obs = adata.obs.merge(
        dcp_tsv,
        left_on=cxg_col,
        right_on=id_col,
        how='left'
    )

    # make unique per barcode
    adata.obs = adata.obs.drop_duplicates(subset=["index"])
    del adata.obs['index']
    assert adata.n_obs == adata.obs.shape[0], f'Number of observations changed from {adata.n_obs} to {adata.obs.shape[0]}'

# save stats
intersect_max.to_csv(out_stats, sep='\t', index=True)
# TODO: metadata column completeness
# TODO: aggregatedness of donor ID

# save anndata
logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obs'],
)
