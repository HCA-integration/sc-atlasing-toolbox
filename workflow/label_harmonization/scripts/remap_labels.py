import pandas as pd
import scanpy as sc

from utils.io import read_anndata

input_anndata = snakemake.input.anndata
input_mapping = snakemake.input.mapping
output_file = snakemake.output.h5ad
mapping_order = snakemake.params.mapping_order


print('mapping order:', mapping_order)

print(f'read adata...')
adata = read_anndata(input_anndata)

label_mapping = pd.read_table(input_mapping)
label_key = None

for mapping_label in mapping_order:

    if label_key is None:
        try:
            assert mapping_label in adata.obs.columns
        except AssertionError:
            raise ValueError(
                f'"{mapping_label}" not found in adata.obs.columns. '
                f'Please make sure the first entry in the mapping order is a column in adata.obs.'
            )
        label_key = mapping_label
        continue

    print(f'mapping "{label_key}" to "{mapping_label}"...')
    df = label_mapping[[mapping_label, label_key]].drop_duplicates()
    map_dict = df.set_index(label_key)[mapping_label].to_dict()
    
    # map...
    mapped = adata.obs[label_key].map(map_dict)
    adata.obs[mapping_label] = pd.Series(mapped, dtype="category")

    # set current mapping label as new label key
    label_key = mapping_label

print('write...')
adata.write_h5ad(output_file)