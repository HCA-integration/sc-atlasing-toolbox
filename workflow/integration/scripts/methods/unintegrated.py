from pathlib import Path
import scanpy as sc

from utils import add_metadata, read_anndata, process, select_layer
from utils_pipeline.io import link_zarr

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

if 'connectivities' not in adata.obsp \
    or 'distances' not in adata.obsp \
    or 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata)
    files_to_keep.extend(['obsp'])

adata = process(adata=adata, adata_raw=adata, output_type=params['output_type'])
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