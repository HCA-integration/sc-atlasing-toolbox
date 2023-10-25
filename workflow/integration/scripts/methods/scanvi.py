import scvi

from utils import add_metadata
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer
from utils_pipeline.processing import process


input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
label_key = wildcards.label

adata_raw = read_anndata(input_adata)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# run method
# adata = scib.ig.scanvi(adata_raw, batch=wildcards.batch, labels=wildcards.label, **params['hyperparams'])

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
train_params = ['max_epochs', 'observed_lib_size', 'n_samples_per_label']
model_params = {k: v for k, v in hyperparams.items() if k not in train_params}
train_params = {k: v for k, v in hyperparams.items() if k in train_params}

# prepare data for model
adata = adata_raw.copy()
adata.layers['counts'] = select_layer(adata, params['raw_counts'])
adata.obs[label_key] = adata.obs[label_key].astype(str).astype('category')

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
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
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)
