import logging
logging.basicConfig(level=logging.INFO)
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
subsample = snakemake.params.get('subsample', False)
force_scale = snakemake.params.get('force_scale', False)
input_layer = snakemake.params.get('input_layer', 'X')

logging.info(f'params: {params}')

logging.info(f'Read {input_file}...')
kwargs = {'obs': 'obs', 'var': 'var', 'uns': 'uns', 'obsm': 'obsm'}
if use_pct:
    kwargs |= {'X': input_layer}
adata = read_anndata(input_file, **kwargs)
print(adata, flush=True)

# subsample genes
if subsample:
    assert 0 < subsample < 1
    sc.pp.subsample(adata, fraction=subsample)

# scale for PCT if not already scaled
scaled = adata.uns.get('preprocessing', {}).get('scaled', False)
if use_pct and (force_scale or not scaled):
    logger.info("Scale .X for PCT")
    adata, _ = subset_hvg(adata, 'highly_variable')
    sc.pp.scale(adata, max_value=10)

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
print(adata.obs, flush=True)
adata.obs[['reannotation_index', 'low_hierarchy', 'high_hierarchy']] = \
    reannotation.loc[adata.obs_names, ['reannotation_index', 'reannotation', 'group']]
print(adata.obs, flush=True)

logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)