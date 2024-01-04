import torch
from scarches.models.scpoli import scPoli
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes
from utils_pipeline.io import read_anndata, write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
label_key = wildcards.label

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
model_params = hyperparams.get('model', {})
train_params = hyperparams.get('train', {})
n_epochs = model_params.get('n_epochs')
pretrain_epochs = int(0.8 * n_epochs) if n_epochs is not None else None
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)


# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/raw_counts',
    obs='obs',
    var='var',
    uns='uns'
)

# prepare data for model
adata.X = adata.X.astype('float32')
if label_key in adata.obs.columns:
    adata.obs[label_key] = adata.obs[label_key].astype(str).fillna('NA').astype('category')

logging.info(f'Set up scPoli with parameters:\n{pformat(model_params)}')
model = scPoli(
    adata=adata,
    condition_keys=[batch_key],
    cell_type_keys=[label_key] if hyperparams.get('supervised', False) else None,
    unknown_ct_names=['NA'],
    **model_params,
)

logging.info(f'Train scPoli with parameters:\n{pformat(train_params)}')
model.train(
    **train_params,
    pretraining_epochs=pretrain_epochs,
    # alpha_epoch_anneal=100,
    early_stopping_kwargs=early_stopping_kwargs,
    batch_size=32,
)

logging.info('Save model...')
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent(adata, mean=True)
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
    train=model.trainer.logs['epoch_loss'],
    validation=model.trainer.logs['val_loss'],
    output_path=f'{output_model}/train_loss.png'
)

logging.info(adata.__str__())
logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['X', 'obsm', 'var', 'varm', 'varp', 'uns']
)
