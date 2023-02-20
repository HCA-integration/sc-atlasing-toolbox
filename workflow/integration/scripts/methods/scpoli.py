from scipy.sparse import csr_matrix
import torch
from scarches.models.scpoli import scPoli

from utils import add_metadata, read_anndata, process

input_file = snakemake.input.h5ad
output_file = snakemake.output.h5ad
output_model = snakemake.output.model
wildcards = snakemake.wildcards
params = snakemake.params

hyperparams = {} if params['hyperparams'] is None else params['hyperparams']
cell_type_keys = [wildcards.label] if 'supervised' in hyperparams.keys() and hyperparams['supervised'] else None
early_stopping_kwargs = hyperparams['early_stopping_kwargs'] if 'early_stopping_kwargs' in hyperparams.keys() else {}

# check GPU
print('GPU available:', torch.cuda.is_available())
# scvi.settings.batch_size = 32

adata_raw = read_anndata(input_file)

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# prepare anndata for training
adata = adata_raw.copy()
adata.X = adata.layers['counts'].todense()

# train model
model = scPoli(
    adata=adata,
    condition_key=wildcards.batch,
    cell_type_keys=cell_type_keys,
    embedding_dim=3,
)

model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)

model.save(output_model, overwrite=True)

# prepare output adata
adata = adata_raw.copy()
adata.obsm["X_emb"] = csr_matrix(model.get_latent())
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write(output_file)
