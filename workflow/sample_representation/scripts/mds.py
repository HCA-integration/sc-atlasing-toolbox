import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked

def distances_to_mds(distances, random_state=42):
    from sklearn.manifold import MDS

    mds = MDS(
        n_components=2,
        random_state=random_state,
        dissimilarity="precomputed"
    )
    coordinates = mds.fit_transform(distances)
    return coordinates


input_file = snakemake.input.zarr
output_file = snakemake.output.zarr

logging.info(f'Read "{input_file}"...')
adata = read_anndata(
    input_file,
    X='obsp/distances',
    obs='obs',
    obsm='obsm',
)

distances = (adata.X + adata.X.T) / 2
adata.obsm['X_mds'] = distances_to_mds(distances.toarray())

logging.info(f'Write "{output_file}"...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm'],
)
