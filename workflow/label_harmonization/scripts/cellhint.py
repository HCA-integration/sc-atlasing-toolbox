import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad
import scanpy as sc
import cellhint

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_reannotation = snakemake.output.reannotation
output_relation = snakemake.output.relation
output_model = snakemake.output.model

author_label_key = snakemake.params.author_label_key
dataset_key = snakemake.params.dataset_key
params = snakemake.params.get('params', {})
use_pct = params.get('use_pct', False)
use_rep = params.get('use_rep', 'X_pca')
subsample = snakemake.params.get('subsample', False)
force_scale = snakemake.params.get('force_scale', False)
input_layer = snakemake.params.get('input_layer', 'X')

if use_rep is None:
    use_rep = 'X_pca'
    params['use_rep'] = use_rep

logging.info(f'params: {params}')

logging.info(f'Read {input_file}...')
kwargs = {'obs': 'obs', 'var': 'var', 'uns': 'uns', 'obsm': 'obsm'}
if use_pct or not input_file.endswith(('.zarr', '.zarr/')):
    kwargs |= {'X': input_layer}
adata = read_anndata(
    input_file,
    **kwargs,
    dask=True,
    backed=True,
)
print(adata, flush=True)
obs = adata.obs

# subsample cells
if subsample:
    assert 0 < subsample < 1
    logging.info(f'Subsample to {subsample} cells...')
    # sc.pp.subsample(adata, fraction=subsample)
    subset_indices = np.random.choice(
        adata.obs_names,
        size=int(adata.n_obs * subsample),
        replace=False
    )
    adata = adata[adata.obs_names.isin(subset_indices)].copy()

subsetted = False

# scale for PCT if not already scaled
scaled = adata.uns.get('preprocessing', {}).get('scaled', False)
if use_pct and (force_scale or not scaled):
    adata, subsetted = subset_hvg(adata, var_column='highly_variable')
    logger.info("Scale .X for PCT")
    sc.pp.scale(adata, max_value=10)

if use_rep == 'X_pca' and use_rep not in adata.obsm.keys():
    logging.info('Compute PCA...')
    if not subsetted:
        adata, _ = subset_hvg(adata, var_column='highly_variable')
    sc.pp.pca(adata, use_highly_variable=True)

logging.info('Harmonize with CellHint...')
alignment = cellhint.harmonize(
    adata=adata,
    dataset=dataset_key,
    cell_type=author_label_key,
    reannotate=True,
    **params
)
alignment.write(output_model)

# Add index to reannotations
reannotation = alignment.reannotation
mapping = {k: str(i) for i, k in enumerate(reannotation['reannotation'].unique())}
reannotation['reannotation_index'] = reannotation['reannotation'].map(mapping)
reannotation.to_csv(output_reannotation, sep='\t')

# Save relations
relation = alignment.relation
relation['group'] = alignment.groups
relation.to_csv(output_relation, sep='\t', index=False)

# add reannotations to adata.obs
anno_cols = ['reannotation_index', 'reannotation', 'group']
if input_file.endswith(('.zarr', '.zarr/')):
    obs[anno_cols] = reannotation.loc[adata.obs_names, anno_cols]
    adata = ad.AnnData(obs=obs)
else:
    adata.obs[anno_cols] = reannotation.loc[adata.obs_names, anno_cols]
print(adata.obs, flush=True)

logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)