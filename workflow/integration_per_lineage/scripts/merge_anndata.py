from pathlib import Path
import mudata as mu
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, link_zarr


input_files = snakemake.input
output_file = snakemake.output[0]
lineages = snakemake.input.keys()

logging.info('Reading data...')
adatas = [read_anndata(file) for file in input_files]
for ad in adatas:
    if 'full' in ad.uns['integration']['output_type']:
        assert not isinstance(ad.X, type(None))

# sort genes
adatas = [ad[:, ad.var.index.sort_values()] for ad in adatas]

mdata = mu.MuData(
    {str(lineage): adata for lineage, adata in zip(lineages, adatas)}
)
mdata.uns = adatas[0].uns

# write file
mdata.write_zarr(output_file)
for lineage, input_file in zip(lineages, input_files):
    if input_file.endswith('.zarr'):
        input_zarr_files = [f.name for f in Path(input_file).iterdir()]
        files_to_link = [f for f in input_zarr_files if f not in ['.zattrs', '.zgroup', 'var']]
        link_zarr(
            in_dir=input_file,
            out_dir=Path(output_file) / 'mod' / lineage,
            file_names=files_to_link,
            overwrite=True,
        )
