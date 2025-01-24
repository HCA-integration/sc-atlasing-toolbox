import torch
import drvi
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    DRVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
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

torch.set_float32_matmul_precision('medium')
# drvi.settings.seed = params.get('seed', 0)
# drvi.settings.progress_bar_style = 'tqdm'
# drvi.settings.num_threads = snakemake.threads

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=DRVI_MODEL_PARAMS,
)
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

# check GPU
print(f'GPU available: {torch.cuda.is_available()}', flush=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/counts',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

# subset features
adata, subsetted = subset_hvg(adata, var_column='integration_features')

if isinstance(categorical_covariate_keys, list):
    for cov in categorical_covariate_keys:
        adata.obs[cov] = adata.obs[cov].astype('str')
else:
    categorical_covariate_keys = []

# add batch key to categorical covariates
if batch_key not in categorical_covariate_keys:
    categorical_covariate_keys.append(batch_key)

drvi.model.DRVI.setup_anndata(
    adata,
    categorical_covariate_keys=categorical_covariate_keys,
    continuous_covariate_keys=continuous_covariate_keys,
    is_count_data=False,
)

model_params |= dict(categorical_covariates=categorical_covariate_keys)
logging.info(f'Set up DRVI with parameters:\n{pformat(model_params)}')
model = drvi.model.DRVI(adata, **model_params)

# train_params |= dict(plan_kwargs={
#     "n_epochs_kl_warmup": train_params.get('max_epochs', 400),
# })
logging.info(f'Train DRVI with parameters:\n{pformat(train_params)}')
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
