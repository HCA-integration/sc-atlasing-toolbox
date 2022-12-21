from pathlib import Path
import scanpy as sc


input_file = snakemake.input.h5ad
output_dir = snakemake.output[0]
split_key = snakemake.params.split_key

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

adata = sc.read(input_file)

for split in adata.obs[split_key]:
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"{split_file}.h5ad"
    adata[adata.obs[split_key] == split].write(out_file)
