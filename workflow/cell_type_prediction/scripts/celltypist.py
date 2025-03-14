import matplotlib.pyplot as plt
import scanpy as sc
import celltypist
from celltypist import models

from utils.io import read_anndata

sc.settings.verbosity = 3
sc.set_figure_params(frameon=False)


input_file = snakemake.input.h5ad
input_model = snakemake.input.model
method = snakemake.wildcards.method
model_name = snakemake.wildcards.model
label_key = snakemake.params.label_key
output = snakemake.output.tsv
output_png = snakemake.output.png

print(f'Read file: {input_file}...', flush=True)
adata = read_anndata(input_file, X='X', obs='obs', var='var')

# assign and subset adata based on feature names instead on ensembl IDs
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']
non_ENS_genes_list = [name for name in adata.var_names if not name.startswith('ENS')]
adata = adata[:, non_ENS_genes_list].copy()

# normalise and log1p the raw counts TODO: move to preprocessing module
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# run celltypist
# Select the model from the above list. If the `model` argument is not provided, will default to `Immune_All_Low.pkl`.
model = models.Model.load(model=input_model)
# The model summary information.
print(model, flush=True)

# Predict the identity of each input cell with the new model.
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
print(predictions, flush=True)

# plot predictions vs author labels
plt.rcParams['figure.figsize'] = 10, 20  # TODO: make dependent on number of cell types
fig, axes = plt.subplots(nrows=2, ncols=1)
celltypist.dotplot(
    predictions,
    use_as_reference=label_key,
    use_as_prediction='majority_voting',
    title=f'CellTypist label transfer {label_key} vs. majority_voting',
    show=False,
    ax=axes[0]
)
celltypist.dotplot(
    predictions,
    use_as_reference=label_key,
    use_as_prediction='predicted_labels',
    title=f'CellTypist label transfer {label_key} vs. predicted_labels',
    show=False,
    ax=axes[1]
)
plt.tight_layout()
plt.savefig(output_png)

# Get an `AnnData` with predicted labels and confidence scores embedded into the observation metadata columns.
prefix = f'{method}_{model_name}:'
results = predictions.to_adata(insert_labels=True, insert_conf=True, prefix=prefix).obs
results = results[[x for x in results.columns if x.startswith(prefix)]]

results.to_csv(output, sep='\t')
