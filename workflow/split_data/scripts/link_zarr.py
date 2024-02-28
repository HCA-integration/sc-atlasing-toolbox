from pathlib import Path
from utils.io import link_zarr

input_dir = Path(snakemake.input[0])
output_file = Path(snakemake.output[0])
wildcards = snakemake.wildcards

if not output_file.exists():
    output_file.mkdir()

input_file = input_dir / f'value~{wildcards.value}.zarr'
assert input_file.exists(), f'File {input_file} does not exist'

link_zarr(
    in_dir=input_file,
    out_dir=output_file,
    overwrite=True,
)