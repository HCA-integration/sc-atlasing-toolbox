import pandas as pd
import anndata
import scanpy as sc

input_anndata = snakemake.input.anndata
input_mapping = snakemake.input.mapping
output_file = snakemake.output.zarr
mapping_order = snakemake.params.mapping_order


def read_anndata(file):
    if file.endswith('.zarr'):
        adata = anndata.read_zarr(file)
    else:
        adata = sc.read(file)
    return adata


adata = read_anndata(input_anndata)
label_mapping = pd.read_table(input_mapping)
label_key = None

print('mapping order:', mapping_order)

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

    adata.obs[mapping_label] = adata.obs[label_key].map(
        df.set_index(label_key)[mapping_label]
    )

    label_key = mapping_label

adata.write_zarr(output_file)