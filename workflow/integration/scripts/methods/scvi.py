import scvi
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots, set_model_history_dtypes
from utils_pipeline.io import read_anndata, write_zarr_linked
from utils_pipeline.accessors import select_layer, subset_hvg


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers', raw='raw')
adata.X = select_layer(adata, params['raw_counts']) # TODO select layer before reading

# subset to HVGs
adata = subset_hvg(adata)

# run method
hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
train_params = ['max_epochs', 'observed_lib_size', 'n_samples_per_label']
model_params = {k: v for k, v in hyperparams.items() if k not in train_params}
train_params = {k: v for k, v in hyperparams.items() if k in train_params}

scvi.model.SCVI.setup_anndata(
    adata,
    # categorical_covariate_keys=[batch], # Does not work with scArches
    batch_key=wildcards.batch,
)

model = scvi.model.SCVI(
    adata,
    **model_params
)
model.train(**train_params)
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation()
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    # history is not saved with standard model saving
    model_history=set_model_history_dtypes(model.history)
)

logging.info(adata.__str__())
logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['X', 'obsm', 'uns']  # TODO: link to correct .X slot?
)