import scvi
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    SCVI_MODEL_PARAMS, plot_model_history
from utils.io import read_anndata, write_zarr_linked, to_memory
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
var_mask = params.get('var_mask', 'highly_variable')

scvi.settings.seed = params.get('seed', 0)
scvi.settings.progress_bar_style = 'tqdm'
scvi.settings.num_threads = snakemake.threads

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCVI_MODEL_PARAMS,
)
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/raw_counts',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)

# prepare data for model
adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')

# subset features
adata, _ = subset_hvg(adata, var_column=var_mask)

scvi.model.SCVI.setup_anndata(
    adata,
    # categorical_covariate_keys=[batch], # Does not work with scArches
    batch_key=batch_key,
)

logging.info(f'Set up scVI with parameters:\n{pformat(model_params)}')
model = scvi.model.SCVI(
    adata,
    **model_params
)

logging.info(f'Train scVI with parameters:\n{pformat(train_params)}')
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

for loss in ['reconstruction_loss', 'elbo', 'kl_local']:
    train_key = f'{loss}_train'
    validation_key = f'{loss}_validation'
    if train_key not in model.history or validation_key not in model.history:
        continue
    plot_model_history(
        title=loss,
        train=model.history[train_key][train_key],
        validation=model.history[validation_key][validation_key],
        output_path=f'{output_plot_dir}/{loss}.png'
    )


logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
