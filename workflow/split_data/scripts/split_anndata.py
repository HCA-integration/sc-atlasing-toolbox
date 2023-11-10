import sys
from pathlib import Path
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import ensure_sparse

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    obsm='obsm',
    varm='varm',
    layers='layers'
)

for layer in ['X']+list(adata.layers.keys()):
    ensure_sparse(adata, layer=layer)

splits = adata.obs[split_key].astype(str).unique()
logging.info(f'splits: {splits}')

for split in splits:
    logging.info(f'Split by {split_key}={split}')
    # split anndata
    adata_sub = adata[adata.obs[split_key] == split]
    # write to file
    split_file = split.replace(' ', '_').replace('/', '_')
    out_file = out_dir / f"value~{split_file}.zarr"

    logging.info(f'write to {out_file}...')
    adata_sub.write_zarr(out_file)
    del adata_sub
