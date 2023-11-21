from pathlib import Path
import pandas as pd

from utils.io import read_anndata, link_zarr


input_anndata = snakemake.input[0]
input_scrublet = snakemake.input.scrublet
input_doubletdetection = snakemake.input.doubletdetection
output_zarr = snakemake.output.zarr


# read AnnData
adata = read_anndata(input_anndata, obs='obs')

scrub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_scrublet])
doub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_doubletdetection])

print(scrub_scores)
print(doub_scores)

adata.obs = adata.obs.join(scrub_scores)
print(adata.obs)
adata.obs = adata.obs.join(doub_scores)
print(adata.obs)

adata.write_zarr(output_zarr)

if input_anndata.endswith('.zarr'):
    input_files = [f.name for f in Path(input_anndata).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obs']]
    link_zarr(
        in_dir=input_anndata,
        out_dir=output_zarr,
        file_names=files_to_keep,
        overwrite=True,
    )