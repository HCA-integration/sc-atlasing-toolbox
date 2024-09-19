import os
from celltypist import models

model = snakemake.wildcards.model
output_model = snakemake.output.model

model_name = f'{model}.pkl'
print(f'model name: {model_name}')
models.download_models(model=model_name, force_update=True)
download_path = f'{snakemake.params["CELLTYPIST_FOLDER"]}/data/models/{model_name}'
os.link(download_path, output_model)
