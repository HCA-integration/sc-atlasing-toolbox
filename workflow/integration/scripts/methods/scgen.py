from scgen import SCGEN

from utils import add_metadata, read_anndata, process

input_file = snakemake.input.h5ad
output_file = snakemake.output.h5ad
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
cell_type_keys = [wildcards.label] if 'supervised' in hyperparams.keys() and hyperparams['supervised'] else None
early_stopping_kwargs = hyperparams['early_stopping_kwargs'] if 'early_stopping_kwargs' in hyperparams.keys() else {}

adata_raw = read_anndata(input_file)
adata = adata_raw.copy()
adata.X = adata.layers['counts']

SCGEN.setup_anndata(adata, batch_key=wildcards.batch, labels_key=wildcards.label)

# train model
model = SCGEN(adata)
model.train()
corrected_adata = model.batch_removal()
model.save(output_model, overwrite=True)

# prepare output adata
adata.X = corrected_adata.X
adata.obsm["X_emb"] = corrected_adata.obsm["latent"]
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_file)
