from scipy.sparse import issparse
import torch
from scarches.models.scpoli import scPoli

from utils import add_metadata, read_anndata, process, select_layer

input_file = snakemake.input.h5ad
output_file = snakemake.output.h5ad
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
cell_type_keys = [wildcards.label] if 'supervised' in hyperparams.keys() and hyperparams['supervised'] else None
model_params = hyperparams['model'] if 'model' in hyperparams.keys() else {}
train_params = hyperparams['train'] if 'train' in hyperparams.keys() else {}
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

adata_raw = read_anndata(input_file)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# prepare anndata for training
adata = adata_raw.copy()
adata.X = select_layer(adata, params['raw_counts'], force_dense=True)
adata.X = adata.X.astype('float32')

# train model
model = scPoli(
    adata=adata,
    condition_keys=wildcards.batch,
    cell_type_keys=cell_type_keys,
    **model_params,
)

model.train(
    **train_params,
    pretraining_epochs=4,
    alpha_epoch_anneal=100,
    early_stopping_kwargs=early_stopping_kwargs,
)

model.save(output_model, overwrite=True)

# prepare output adata
adata = adata_raw.copy()
adata.obsm["X_emb"] = model.get_latent(adata, mean=True)
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_file, compression='lzf')
