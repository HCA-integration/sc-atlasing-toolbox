import sys
from pathlib import Path
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad

from utils.io import read_anndata
from utils.misc import ensure_sparse

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key
values = snakemake.params.get('values', [])

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file, backed=True)

for layer in ['X']+list(adata.layers.keys()):
    ensure_sparse(adata, layer=layer)

file_value_map = {
    s.replace(' ', '_').replace('/', '_'): s
    for s in adata.obs[split_key].astype(str).unique()
}
split_files = list(file_value_map.keys())
split_files = set(split_files + values)
logging.info(f'splits: {split_files}')

for split_file in split_files:
    split = file_value_map.get(split_file)
    logging.info(f'Split by {split_key}={split}')
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split]
    # write to file
    out_file = out_dir / f"value~{split_file}.zarr"

    logging.info(f'write to {out_file}...')
    if adata_sub.n_obs == 0:
        adata_sub = ad.AnnData(obs=adata_sub.obs, var=adata_sub.var)
    adata_sub.write_zarr(out_file)
    del adata_sub
