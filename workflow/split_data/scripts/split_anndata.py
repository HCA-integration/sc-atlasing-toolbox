import faulthandler
faulthandler.enable()
from pathlib import Path
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad
from pprint import pformat
from scipy.sparse import csr_matrix, coo_matrix
import sparse
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, csr_matrix_int64_indptr
from utils.accessors import adata_to_memory
from utils.misc import apply_layers

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key
values = snakemake.params.get('values', [])
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)
exclude_slots = snakemake.params.get('exclude_slots', [])

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

# convert split_key column to string
adata.obs[split_key] = adata.obs[split_key].astype(str)
logging.info(adata.obs[split_key].value_counts())

# logging.info('Convert dtypes...')
# def convert_to_dtype(x, func):
#     if isinstance(x, da.Array):
#         return x.map_blocks(lambda x: func(x), dtype=x.dtype)
#     # if isinstance(x, csr_matrix):
#     #     return func(x)
#     return x

# adata = apply_layers(
#     adata,
#     lambda x: convert_to_dtype(x, sparse.GCXS)
# )

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
        adata_sub = adata_sub.copy()
        
        if backed:
            logging.info('Load backed arrays to memory...')
            adata_sub = adata_to_memory(adata_sub)
        if dask:
            logging.info('Convert to csr_matrix with int64 indptr...')
            adata_sub = apply_layers(
                adata_sub,
                # lambda x: convert_to_dtype(x, lambda y: y.to_scipy_sparse())
                # lambda x: x.compute().to_scipy_sparse() if isinstance(x, da.Array) else x
                # lambda x: csr_matrix(x.compute()) if isinstance(x, da.Array) else x
                lambda x: x.map_blocks(csr_matrix_int64_indptr, dtype=x.dtype)
            )
    
    # write to file
    logging.info(f'Write to {out_file}...')
    adata_sub.write_zarr(out_file)
    del adata_sub
