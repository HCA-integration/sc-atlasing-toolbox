import sys
from pathlib import Path
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad
from pprint import pformat
from dask import array as da

from utils.io import read_anndata
from utils.accessors import adata_to_memory
from utils.misc import apply_layers
from utils.sparse_dask import sparse_dataset_as_dask

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key
values = snakemake.params.get('values', [])
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)
exclude_slots = snakemake.params.get('exclude_slots', [])
n_threads = snakemake.threads

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(
    input_file,
    backed=backed,
    dask=dask,
    exclude_slots=exclude_slots,
)
logging.info(adata.__str__())

file_value_map = {
    s.replace(' ', '_').replace('/', '_'): s
    for s in adata.obs[split_key].astype(str).unique()
}
logging.info(f'file_value_map: {pformat(file_value_map)}')
# split_files = list(file_value_map.keys())
# split_files = set(split_files + values)
split_files = values
logging.info(f'splits: {split_files}')

for split_file in split_files:
    split = file_value_map.get(split_file)
    out_file = out_dir / f"value~{split_file}.zarr"
    
    # split anndata
    logging.info(f'Split by {split_key}={split}')
    adata_sub = adata[adata.obs[split_key] == split]
    logging.info(adata_sub.__str__())
    
    if adata_sub.n_obs == 0:
        adata_sub = ad.AnnData(
            X=np.zeros(adata_sub.shape),
            obs=adata_sub.obs,
            obsm=adata_sub.obsm,
            obsp=adata_sub.obsp,
            var=adata_sub.var,
            varm=adata_sub.varm,
            varp=adata_sub.varp,
            layers={
                layer: np.zeros(adata_sub.shape)
                for layer in adata_sub.layers
            },
            uns=adata_sub.uns,
        )
    else:
        logging.info('Copy subset...')
        adata_sub = apply_layers(
            adata_sub.copy(),
            lambda x: x.rechunk({0: 1000}) if isinstance(x, da.Array) else x
        )
        adata_sub = adata_to_memory(adata_sub)
    
    # write to file
    logging.info(f'Write to {out_file}...')
    adata_sub.write_zarr(out_file)
    del adata_sub
