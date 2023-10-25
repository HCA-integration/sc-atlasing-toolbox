from pathlib import Path
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr
from utils_pipeline.accessors import select_layer


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata = read_anndata(input_file)
adata.X = select_layer(adata, params['norm_counts'])

# prepare output adata
files_to_keep = ['obsm', 'uns', 'layers']

if 'X_pca' not in adata.obsm:
    sc.pp.pca(adata, use_highly_variable=True)
    files_to_keep.extend(['varm'])
adata.obsm['X_emb'] = adata.obsm['X_pca']

logging.info(adata.__str__())
logging.info(adata.uns)
if 'connectivities' not in adata.obsp \
    or 'distances' not in adata.obsp \
    or 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata)
else:
    logging.info(adata.uns['neighbors'])
    logging.info(adata.obsp.keys())
    files_to_keep.extend(['obsp', 'uns'])

adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# write file
del adata.layers
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_link = [f for f in input_files if f not in files_to_keep]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_link,
        overwrite=True,
    )
