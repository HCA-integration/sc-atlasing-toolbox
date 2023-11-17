from pathlib import Path
from utils.io import link_zarr

input_file = Path(snakemake.input[0])
output_file = Path(snakemake.output[0])

if not output_file.exists():
    output_file.mkdir()

link_zarr(
    in_dir=input_file,
    out_dir=output_file,
    overwrite=True,
)