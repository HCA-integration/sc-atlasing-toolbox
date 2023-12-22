import scvi
import logging
import pickle
import os
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial
from utils_pipeline.accessors import select_layer

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
train_params = ['max_epochs', 'observed_lib_size', 'n_samples_per_label', 'batch_size']
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
model.train(**train_params, early_stopping=True)
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation()
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(
    adata,
    wildcards,
    params,
    # history is not saved with standard model saving
    model_history=set_model_history_dtypes(model.history)
)
