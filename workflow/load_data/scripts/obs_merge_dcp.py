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

obs_df = read_elem(zarr.open(in_file)['obs'])
dcp_tsv = pd.read_table(in_dcp)

obs_ids = set(obs_df['donor_id'].unique())
n_donors = len(obs_ids)

# identify ID column
cols = []
intersection = []

for col in id_cols:
    if col not in dcp_tsv.columns:
        continue
    
    # explode ID column
    dcp_tsv_exploded = explode_table(dcp_tsv, col)

    # count ID overlaps
    dcp_ids = set(dcp_tsv_exploded[col].unique())
    intersect = obs_ids.intersection(dcp_ids)
    n_intersect = len(intersect)

    cols.append(col)
    intersection.append(n_intersect)

intersect_df = pd.DataFrame({'columns': cols, 'intersection': intersection}).set_index('columns')
intersect_df['cxg_donors'] = n_donors
intersect_df['intersection_fraction'] = intersect_df['intersection'] / n_donors
print(intersect_df)

argmax = intersect_df['intersection'].argmax()
id_col = intersect_df.index.tolist()[argmax]

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
    left_on='donor_id',
    right_on=id_col,
    how='left'
)

# make unique per barcode
obs_df = obs_df.drop_duplicates(subset=["barcode"])

# save obs
print(obs_df)
obs_df.to_csv(out_obs, sep='\t', index=False)

# save stats
intersect_df.loc[id_col,:][['intersection', 'cxg_donors', 'intersection_fraction']].to_csv(out_stats, sep='\t', index=True)
# TODO: metadata column completeness
# TODO: aggregatedness of donor ID
