"""
Copied and adapted from https://gist.github.com/ivirshup/c29c9fb0b5b21a9c290cf621e4e68b18
"""
import h5py
from scipy import sparse
import dask.array as da
from dask import delayed
import zarr
import anndata as ad
from anndata.experimental import read_elem, sparse_dataset


def csr_callable(shape: tuple[int, int], dtype) -> sparse.csr_matrix:
    if len(shape) == 0:
        shape = (0, 0)
    if len(shape) == 1:
        shape = (shape[0], 0)
    elif len(shape) == 2:
        pass
    else:
        raise ValueError(shape)
    return sparse.csr_matrix(shape, dtype=dtype)


class CSRCallable:
    """Dummy class to bypass dask checks"""
    def __new__(cls, shape, dtype):
        return csr_callable(shape, dtype)


def make_dask_chunk(x: "SparseDataset", start: int, end: int) -> da.Array:
    def take_slice(x, idx):
        return x[idx]

    return da.from_delayed(
        delayed(take_slice)(x, slice(start, end)),
        dtype=x.dtype,
        shape=(end - start, x.shape[1]),
        meta=CSRCallable,
    )


def sparse_dataset_as_dask(x, stride: int):
    n_chunks, rem = divmod(x.shape[0], stride)

    chunks = []
    cur_pos = 0
    for i in range(n_chunks):
        chunks.append(make_dask_chunk(x, cur_pos, cur_pos + stride))
        cur_pos += stride
    if rem:
        chunks.append(make_dask_chunk(x, cur_pos, x.shape[0]))

    return da.concatenate(chunks, axis=0)


def read_w_sparse_dask(group: [h5py.Group, zarr.Group], obs_chunk: int = 1000) -> ad.AnnData:
    return ad.AnnData(
        X=sparse_dataset_as_dask(sparse_dataset(group["X"]), obs_chunk),
        **{
            k: read_elem(group[k]) if k in group else {}
            for k in ["layers", "obs", "var", "obsm", "varm", "uns", "obsp", "varp"]
        }
    )
