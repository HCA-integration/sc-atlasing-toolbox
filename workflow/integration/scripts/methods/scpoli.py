from scipy.sparse import issparse
import torch
from scarches.models.scpoli import scPoli

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial
from utils_pipeline.accessors import select_layer


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
pretrain_epochs = int(0.8 * model_params['n_epochs']) if 'n_epochs' in model_params else None
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}


# check GPU
print('GPU available:', torch.cuda.is_available())
# scvi.settings.batch_size = 32

adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers')
adata.X = select_layer(adata, params['norm_counts'])

# subset to HVGs
adata = adata[:, adata.var['highly_variable']]

# prepare anndata for training
adata.X = select_layer(adata, params['raw_counts'], force_sparse=True)
adata.X = adata.X.astype('float32')

if label_key in adata.obs.columns:
    adata.obs[label_key] = adata.obs[label_key].astype(str).fillna('NA').astype('category')

# train model
model = scPoli(
    adata=adata,
    condition_keys=[batch_key],
    cell_type_keys=[label_key] if hyperparams.get('supervised', False) else None,
    unknown_ct_names=['NA'],
    **model_params,
)

model.train(
    **train_params,
    pretraining_epochs=pretrain_epochs,
    # alpha_epoch_anneal=100,
    early_stopping_kwargs=early_stopping_kwargs,
    batch_size=32,
)

model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent(adata, mean=True)
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])