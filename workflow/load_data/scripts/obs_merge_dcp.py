from pprint import pprint
import pandas as pd
from anndata.experimental import read_elem
import zarr

in_file = snakemake.input[0]
in_dcp = snakemake.input[1]
out_obs = snakemake.output.obs
out_stats = snakemake.output.stats

def explode_table(df, col, sep=' \|\| '):
    df = df.copy()
    df[col] = df[col].astype(str).str.split(pat=sep)
    return df.explode(col)


id_cols = snakemake.params.id_cols

dcp_tsv = pd.read_table(in_dcp)
with zarr.open(in_file) as z:
    obs_df = read_elem(z['obs'])

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
print(intersect_df)

# identify ID column with most overlap
argmax = intersect_df['intersection'].argmax()
intersect_max = intersect_df.iloc[argmax,:]
assert isinstance(intersect_max, pd.Series), f'max intersection is not a Series\n{intersect_max}'
id_col = intersect_max['dcp_column']
cxg_col = intersect_max['cxg_column']

metadata_columns  = [id_col] + snakemake.params.metadata_cols
metadata_columns = [c for c in metadata_columns if c in dcp_tsv.columns]
print('DCP metadata columns:')
pprint(metadata_columns)

# merge on ID column
dcp_tsv = explode_table(dcp_tsv, id_col).dropna(subset=[id_col])
dcp_tsv = dcp_tsv[metadata_columns].drop_duplicates()
print(dcp_tsv)
obs_df = obs_df.merge(
    dcp_tsv,
    left_on=cxg_col,
    right_on=id_col,
    how='left'
)

# make unique per barcode
obs_df = obs_df.drop_duplicates(subset=["barcode"])

# save obs
print(obs_df)
obs_df.to_csv(out_obs, sep='\t', index=False)

# save stats
intersect_max[['dcp_column', 'cxg_column', 'intersection', 'n_cxg', 'intersection_fraction']].to_csv(out_stats, sep='\t', index=True)
# TODO: metadata column completeness
# TODO: aggregatedness of donor ID
