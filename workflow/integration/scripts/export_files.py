from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored")

from utils.io import read_anndata, link_zarr


input_file = snakemake.input[0]
input_clusters = snakemake.input.clusters
output_file = snakemake.output[0]

adata = read_anndata(input_file)
cluster_df = pd.read_table(input_clusters, index_col=0, dtype=str)

# add cluster assignment to adata
cluster_cols = cluster_df.columns
adata.obs[cluster_cols] = cluster_df.loc[adata.obs_names, cluster_cols]

# write file
adata.write_zarr(output_file)

input_files = [f.name for f in Path(input_file).iterdir()]
files_to_keep = [f for f in input_files if f not in ['obs']]
link_zarr(
    in_dir=input_file,
    out_dir=output_file,
    file_names=files_to_keep,
    overwrite=True,
)
