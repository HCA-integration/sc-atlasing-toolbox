import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from pathlib import Path
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity

from utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes
from utils_pipeline.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
batch_key = wildcards.batch
label_key = wildcards.label
params = snakemake.params

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    train_params=[
        'n_epochs',
        'max_epochs',
        'observed_lib_size',
        'n_samples_per_label'
    ],
)
early_stopping_kwargs = {
    "early_stopping_metric": "val_loss",
    "patience": 20,
    "threshold": 0,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
train_params['early_stopping_kwargs'] = early_stopping_kwargs
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/norm_counts',
    obs='obs',
    var='var',
    uns='uns'
)

# prepare data for model
adata.X = adata.X.astype('float32')
adata = remove_sparsity(adata)

logging.info(f'Set up scGEN with parameters:\n{pformat(model_params)}')
model = sca.models.scgen(
    adata=adata,
    **model_params,
)

logging.info(f'Train scGEN with parameters:\n{pformat(train_params)}')
model.train(**train_params)

logging.info('Save model...')
model.save(output_model, overwrite=True)

logging.info(adata.__str__())
corrected_adata = model.batch_removal(
    adata,
    batch_key=batch_key,
    cell_label_key=label_key,
    return_latent=True,
)

# prepare output adata
adata.X = corrected_adata.X
adata.obsm["X_emb"] = corrected_adata.obsm["latent_corrected"]
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    model_history=dict(model.trainer.logs)
)

# plot model history
from utils import plot_model_history

plot_model_history(
    title='loss',
    train=model.trainer.logs['epoch_loss'],
    validation=model.trainer.logs['val_loss'],
    output_path=f'{output_plot_dir}/loss.png'
)

logging.info(adata.__str__())
logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['X', 'obsm', 'var', 'varm', 'varp', 'uns']
)
