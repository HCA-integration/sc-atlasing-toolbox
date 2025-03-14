from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import celltypist
from celltypist import models

from utils.io import read_anndata, write_zarr_linked

sc.settings.verbosity = 3
sc.set_figure_params(frameon=False)


input_file = snakemake.input[0]
input_model = snakemake.input.model
output_file = snakemake.output[0]

model_name = snakemake.wildcards.celltypist_model
label_key = snakemake.params.label_key
output_png = Path(snakemake.output.png)
output_png.mkdir(parents=True, exist_ok=True)


print(f'Read file: {input_file}...', flush=True)
adata = read_anndata(input_file, X='X', obs='obs', var='var')

# assign and subset adata based on feature names instead on ensembl IDs
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']
non_ENS_genes_list = [name for name in adata.var_names if not name.startswith('ENS')]
adata = adata[:, non_ENS_genes_list].copy()

print('Normalizing and log-transforming data...', flush=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# run celltypist
model = models.Model.load(model=input_model)
print(model, flush=True)

# Predict the identity of each input cell with the new model.
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
print(predictions, flush=True)

# plot predictions vs author labels
celltypist.dotplot(
    predictions,
    use_as_reference=label_key,
    use_as_prediction='predicted_labels',
    title=f'CellTypist label transfer {label_key} vs. predicted_labels',
    # show=False,
)
plt.savefig(output_png / 'predicted_labels.png', bbox_inches='tight')

celltypist.dotplot(
    predictions,
    use_as_reference=label_key,
    use_as_prediction='majority_voting',
    title=f'CellTypist label transfer {label_key} vs. majority_voting',
    # show=False,
)
plt.savefig(output_png / 'majority_voting.png', bbox_inches='tight')

# Get an `AnnData` with predicted labels and confidence scores embedded into the observation metadata columns.
prefix = f'celltypist_{model_name}:'
obs = predictions.to_adata(insert_labels=True, insert_conf=True, prefix=prefix).obs
adata.obs = obs[[x for x in obs.columns if x.startswith(prefix)]]

print(f'Write file: {output_file}...', flush=True)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)
