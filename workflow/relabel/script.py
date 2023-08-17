import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import anndata
from anndata.experimental import read_elem
import zarr

from utils.io import link_zarr

input_anndata = snakemake.input.anndata
input_mapping = snakemake.input.mapping
output_file = snakemake.output.zarr
mapping_order = snakemake.params.mapping_order


logging.info(f'Mapping order:\n{mapping_order}')
label_mapping = pd.read_table(input_mapping, comment='#')
label_key = None

logging.info('Read adata...')
if input_anndata.endswith('.h5ad'):
    adata = anndata.read_h5ad(input_anndata, backed='r')
    obs = adata.obs
    link_output = False
elif input_anndata.endswith('.zarr'):
    z = zarr.open(input_anndata)
    obs = read_elem(z["obs"])
    adata = anndata.AnnData(obs=obs)
    link_output = True
else:
    raise ValueError(f'Unknown file extension: "{input_anndata}"')

for mapping_label in mapping_order:

    if label_key is None:
        try:
            assert mapping_label in obs.columns
        except AssertionError as e:
            raise ValueError(
                f'"{mapping_label}" not found in adata.obs.columns. '
                f'Please make sure the first entry in the mapping order is a column in adata.obs.'
            ) from e
        label_key = mapping_label
        continue

    logging.info(f'mapping "{label_key}" to "{mapping_label}"...')
    df = label_mapping[[mapping_label, label_key]].drop_duplicates()
    map_dict = df.set_index(label_key)[mapping_label].to_dict()

    logging.info('map...')
    mapped = obs[label_key].map(map_dict)
    obs[mapping_label] = pd.Series(mapped, dtype="category")

    # set current mapping label as new label key
    label_key = mapping_label

logging.info('Write...')
adata.obs = obs
adata.write_zarr(output_file)

if link_output:
    input_files = [f.name for f in Path(input_anndata).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obs']]
    link_zarr(
        in_dir=input_anndata,
        out_dir=output_file,
        file_names=files_to_keep,
        overwrite=True,
    )
