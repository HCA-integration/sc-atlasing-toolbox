import scvi

from utils import add_metadata, read_anndata, process, select_layer

input_adata = snakemake.input.h5ad
output_adata = snakemake.output.h5ad
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

adata_raw = read_anndata(input_adata)

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# run method
# adata = scib.ig.scvi(adata_raw, batch=wildcards.batch, **params['hyperparams'])

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
train_params = ['max_epochs', 'observed_lib_size']
model_params = {k: v for k, v in hyperparams.items() if k not in train_params}
train_params = {k: v for k, v in hyperparams.items() if k in train_params}

adata = adata_raw.copy()
adata.layers['counts'] = select_layer(adata, params['raw_counts'])
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    # categorical_covariate_keys=[batch], # Does not work with scArches
    batch_key=wildcards.batch,
)

model = scvi.model.SCVI(
    adata,
    n_latent=20,
    # scArches params
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
    **model_params
)
model.train(**train_params)
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation()
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_adata, compression='gzip')
