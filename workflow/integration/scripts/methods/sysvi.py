from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

import os
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.pyplot import rcParams

import sys
import torch
import pytorch_lightning as pl
#path_to_sysvi = ''
#sys.path.insert(0, path_to_sysvi)
from cross_system_integration.model._xxjointmodel import XXJointModel


from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    SCVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked, to_memory
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

torch.set_float32_matmul_precision('medium')
scvi.settings.seed = params.get('seed', 0)
scvi.settings.progress_bar_style = 'tqdm'
scvi.settings.num_threads = snakemake.threads

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCVI_MODEL_PARAMS,
)
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]
logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

# check GPU
print(f'GPU available: {torch.cuda.is_available()}', flush=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/counts',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

# subset features
adata, subsetted = subset_hvg(adata, var_column='integration_features')

if isinstance(categorical_covariate_keys, list):
    for cov in categorical_covariate_keys:
        adata.obs[cov] = adata.obs[cov].astype('str')


def train_sysvi_model(adata, system_key='Dataset', condition_key='Condition', max_epochs=200, embed_save_path='sysVI/sysVI_dataset_reference_latent.h5ad',
                       model_save_path='sysVI/sysVI_model'):
    """
    Trains an XXJointModel on the provided AnnData object and saves the resulting embedding.

    Parameters:
        adata (AnnData): The annotated data matrix.
        system_key (str): Key in adata.obs to use as system labels.
        condition_key (str): Key in adata.obs to use as condition labels.
        max_epochs (int): Maximum number of training epochs.
        embed_save_path (str): Path to save the embedding AnnData object.

    Returns:
        None
    """
    try:
        logging.info('Setting up adata for XXJointModel')
        adata_copy = adata.copy()
        adata_prepared = XXJointModel.setup_anndata(adata_copy, group_key=None, system_key=system_key, categorical_covariate_keys=[condition_key])

        logging.info('Initialising XXJointModel')
        model = XXJointModel(adata=adata_prepared)
                             # prior='vamp', n_prior_components=5, pseudoinputs_data_init=True,
                             # trainable_priors=True,
                             # encode_pseudoinputs_on_eval_mode=True,
                             # z_dist_metric = 'MSE_standard', n_layers=2,
                             # n_hidden=256)

        logging.info('Starting training')
        model.train(
            max_epochs=100,
            log_every_n_steps=1,
            check_val_every_n_epoch=1,
            val_check_interval=1.0,
            # train_size=0.9,
            # plan_kwargs={
            #     'optimizer': "Adam",
            #     'lr': 0.001,
            #     'reduce_lr_on_plateau': False,
            #     'lr_scheduler_metric': 'loss_train',  # Replace with default value
            #     'lr_patience': 5,  # Replace with default value
            #     'lr_factor': 0.1,  # Replace with default value
            #     'lr_min': 1e-7,  # Replace with default value
            #     'lr_threshold_mode': 'rel',  # Replace with default value
            #     'lr_threshold': 0.1,  # Replace with default value
            #     'log_on_epoch': True,  # Replace with default value
            #     'log_on_step': False,  # Replace with default value
            #     'loss_weights': {
            #         'kl_weight': 1.0,  # Replace with default value
            #         'reconstruction_weight': 1.0,  # Replace with default value
            #         'z_distance_cycle_weight': 5.0, 
            #     },
            # }
        )
        logging.info('Generating embedding')
        embed = model.embed(adata=adata_prepared)
        #adata.obsm['X_sysvi'] = embed
        #adata.uns['output_type'] = 'embed'
        adata.obsm['X_emb'] = embed


    except Exception as e:
        print(f'XXJointModel training failed because {e}')

train_sysvi_model(
    adata=adata,
    system_key='ID_batch_covariate',
    condition_key=batch_key,
    max_epochs=200)

logging.info("sysVI training completed.")





adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    # history is not saved with standard model saving
    #model_history=set_model_history_dtypes(model.history)
)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
