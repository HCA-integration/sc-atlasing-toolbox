import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import anndata
from anndata.experimental import read_elem
import zarr

from utils.io import read_anndata, link_zarr_partial

input_file = snakemake.input.anndata
input_new_cols = snakemake.input.get('new_columns')
input_merge_cols = snakemake.input.get('merge_columns')
output_file = snakemake.output.zarr
file_id = snakemake.wildcards.file_id

logging.info('Read adata...')
adata = read_anndata(input_file, obs='obs')


# merge new columns
if input_new_cols is not None:
    mapping_order = snakemake.params.get('mapping_order')
    logging.info(f'Mapping order:\n{mapping_order}')
    label_mapping = pd.read_table(input_new_cols, comment='#')
    
    label_key = None
    for mapping_label in mapping_order:
        if label_key is None:
            try:
                assert mapping_label in adata.obs.columns
            except AssertionError as e:
                raise ValueError(
                    f'"{mapping_label}" not found in adata.obs.columns. '
                    f'Please make sure the first entry in the mapping order is a column in adata.obs.'
                ) from e
            label_key = mapping_label
            continue
        
        logging.info(f'mapping "{label_key}" to "{mapping_label}"...')
        
        # get unique mapping
        df = label_mapping[[label_key, mapping_label]].drop_duplicates()
        
        # remove trailing whitespaces
        remove_trailing_whitespaces = lambda x: x.str.strip() if hasattr(x, 'str') else x
        adata.obs[label_key] = adata.obs[label_key].apply(remove_trailing_whitespaces)
        df[mapping_label] = df[mapping_label].apply(remove_trailing_whitespaces)
        
        # apply mapping
        map_dict = df.set_index(label_key)[mapping_label].to_dict()
        adata.obs[mapping_label] = pd.Series(adata.obs[label_key].map(map_dict), dtype="category")

        # set current mapping label as new label key
        # label_key = mapping_label

# merge existing columns
if input_merge_cols is not None:
    sep = snakemake.params.get('merge_sep', '-')
    logging.info(f'Merge existing columns with sep="{sep}"...')
    
    merge_config_df = pd.read_table(input_merge_cols, comment='#')
    for col in ['file_id', 'column_name', 'columns']:
        assert col in merge_config_df.columns, f'"{col}" not found in {input_merge_cols}\n{merge_config_df}'
    
    merge_config_df = merge_config_df[merge_config_df['file_id'] == file_id]
    assert not merge_config_df.duplicated().any(), f'Duplicated rows in {input_merge_cols} for {file_id}\n{merge_config_df.duplicated()}'
    
    for _, row in merge_config_df.iterrows():
        col_name = row['column_name']
        cols = row['columns'].split(',')
        logging.info(f'merge {col_name} from {cols}...')
        adata.obs[col_name] = adata.obs[cols].apply(
            lambda x: sep.join(x.values.astype(str)), axis=1
        )

logging.info('Write to {output_file}...')
adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obs'])