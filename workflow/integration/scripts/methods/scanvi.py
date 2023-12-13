import scvi

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

adata = read_anndata(input_file, X='X', obs='obs', var='var', layers='layers', raw='raw', uns='uns')
adata.X = select_layer(adata, params['raw_counts'])

# subset to HVGs
adata = adata[:, adata.var['highly_variable']].copy()

# run method
# adata = scib.ig.scanvi(adata, batch=wildcards.batch, labels=wildcards.label, **params['hyperparams'])

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
train_params = ['max_epochs', 'observed_lib_size', 'n_samples_per_label']
model_params = {k: v for k, v in hyperparams.items() if k not in train_params}
train_params = {k: v for k, v in hyperparams.items() if k in train_params}

# prepare data for model
adata.obs[label_key] = adata.obs[label_key].astype(str).astype('category')

scvi.model.SCVI.setup_anndata(
    adata,
    # categorical_covariate_keys=[batch], # Does not work with scArches
    batch_key=batch_key,
)

# train SCVI
model = scvi.model.SCVI(
    adata,
    **model_params
)
model.train(**train_params)

# train SCANVI on top of SCVI
model = scvi.model.SCANVI.from_scvi_model(
    model,
    labels_key=label_key,
    unlabeled_category='nan',
    **model_params
)
model.train(**train_params)
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation()
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_file)
link_zarr_partial(input_file, output_file, files_to_keep=['obsm', 'uns'])