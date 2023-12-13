import logging
logging.basicConfig(level=logging.INFO)
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, link_zarr_partial
from utils_pipeline.accessors import select_layer

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
wildcards = snakemake.wildcards
batch_key = wildcards.batch
label_key = wildcards.label
params = snakemake.params

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
train_params = ['n_epochs', 'max_epochs', 'observed_lib_size', 'n_samples_per_label']
model_params = {k: v for k, v in hyperparams.items() if k not in train_params}
train_params = {k: v for k, v in hyperparams.items() if k in train_params}
early_stopping_kwargs = {
    "early_stopping_metric": "val_loss",
    "patience": 20,
    "threshold": 0,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
logging.info(hyperparams)

# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers', uns='uns')
adata.X = select_layer(adata, params['norm_counts'])

# subset to HVGs
adata = adata[:, adata.var['highly_variable']]

adata.X = adata.X.astype('float32')
adata = remove_sparsity(adata) # remove sparsity

# scgen.SCGEN.setup_anndata(adata, batch_key=wildcards.batch, labels_key=wildcards.label)
# model = scgen.SCGEN(adata)

model = sca.models.scgen(
    adata=adata,
    **model_params,
)

# train model
model.train(
    **train_params,
    early_stopping_kwargs=early_stopping_kwargs,
)

print(adata)
corrected_adata = model.batch_removal(
    adata,
    batch_key=batch_key,
    cell_label_key=label_key,
    return_latent=True,
)
model.save(output_model, overwrite=True)

# prepare output adata
adata.X = corrected_adata.X
adata.obsm["X_emb"] = corrected_adata.obsm["latent_corrected"]
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# Save output
adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])