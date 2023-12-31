import scvi
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes
from utils_pipeline.io import read_anndata, write_zarr_linked
from utils_pipeline.accessors import select_layer, subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    train_params=[
        'max_epochs',
        'observed_lib_size',
        'n_samples_per_label',
        'batch_size',
        'early_stopping'
    ],
)
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers', raw='raw')
adata.X = select_layer(adata, params['raw_counts']) # TODO select layer before reading

# subset to HVGs
adata = subset_hvg(adata)

scvi.model.SCVI.setup_anndata(
    adata,
    # categorical_covariate_keys=[batch], # Does not work with scArches
    batch_key=wildcards.batch,
)

logging.info(f'Set up scVI with model parameters:\n{pformat(model_params)}')
model = scvi.model.SCVI(
    adata,
    **model_params
)

logging.info(f'Train scVI with training parameters:\n{pformat(train_params)}')
model.train(**train_params)

logging.info('Save model...')
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
    files_to_keep=['X', 'obsm', 'var', 'varm', 'varp', 'uns']  # TODO: link to correct .X slot?
)