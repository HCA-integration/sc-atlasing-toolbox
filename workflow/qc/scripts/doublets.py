from pathlib import Path
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, link_zarr

input_zarr = snakemake.input.zarr
output_obs = snakemake.output.obs
output_zarr = snakemake.output.zarr
batch_key = snakemake.params.get('batch_key')

adata = read_anndata(snakemake.input[0], X='X', obs='obs', var='var')

sc.external.pp.scrublet(
    adata,
    batch_key=batch_key,
    sim_doublet_ratio=2.0,
    expected_doublet_rate=0.05,
    stdev_doublet_rate=0.02,
    synthetic_doublet_umi_subsampling=1.0,
    knn_dist_metric='euclidean',
    normalize_variance=True,
    log_transform=False,
    mean_center=True,
    n_prin_comps=30,
    use_approx_neighbors=True,
    get_doublet_neighbor_parents=False,
    n_neighbors=None,
    threshold=None,
    verbose=True,
    copy=False,
    random_state=0
)

adata.obs.to_csv(output_obs, sep='\t')
adata.write_zarr(output_zarr)

if input_zarr.endswith('.zarr'):
    input_files = [f.name for f in Path(input_zarr).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obs']]
    link_zarr(
        in_dir=input_zarr,
        out_dir=output_zarr,
        file_names=files_to_keep,
        overwrite=True,
    )
