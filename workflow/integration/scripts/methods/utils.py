# from utils_pipeline.misc import ensure_sparse
import numpy as np
from scipy import sparse

SCVI_MODEL_PARAMS = [
    'n_hidden',
    'n_latent',
    'n_layers',
    'dropout_rate',
    'dispersion',
    'gene_likelihood',
    'latent_distribution',
    # VAE parameters
    'n_batch',
    'n_labels',
    'n_hidden',
    'n_latent',
    'n_layers',
    'n_continuous_cov',
    'n_cats_per_cov',
    'dropout_rate',
    'dispersion',
    'log_variational',
    'gene_likelihood',
    'latent_distribution',
    'encode_covariates',
    'deeply_inject_covariates',
    'use_batch_norm',
    'use_layer_norm',
    'use_size_factor_key',
    'use_observed_lib_size',
    'library_log_means',
    'library_log_vars',
    'var_activation',
    'extra_encoder_kwargs',
    'extra_decoder_kwargs',
]

SCANVI_MODEL_PARAMS = SCVI_MODEL_PARAMS + ['linear_classifier']


def ensure_sparse(adata):
    from scipy.sparse import csr_matrix, issparse

    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)


def add_metadata(adata, wildcards, params, **kwargs):
    """
    Add integration metatdata to integratd output
    :param adata:
    :param wildcards:
    :param params:
    :return:
    """
    # TODO: transfer parameters from .uns['preprocessing']

    adata.uns['dataset'] = wildcards.dataset

    if 'methods' in adata.uns:
        adata.uns['methods'].append(wildcards.method)
    else:
        adata.uns['methods'] = [wildcards.method]

    adata.uns['integration'] = {
        'method': wildcards.method,
        'label_key': wildcards.label,
        'batch_key': wildcards.batch,
        'output_type': params['output_type'],
        'hyperparams': params['hyperparams'],
        **kwargs
    }


def remove_slots(adata, output_type, keep_X=False):
    """
    Remove slots that are redundant to integration output
    """
    if isinstance(output_type, str):
        output_type = [output_type]
    
    del adata.layers
    #if 'X_pca' in adata.obsm:
    #    del adata.obsm['X_pca']
    
    if 'full' in output_type:
        ensure_sparse(adata)
    elif 'embed' in output_type:
        if keep_X:
            assert adata.X is not None
            ensure_sparse(adata)
        else:
            del adata.X
    elif 'knn' in output_type:
        # del adata.X
        pass
    else:
        raise ValueError(f'Invalid output type {output_type}')
    return adata


def set_model_history_dtypes(model_history, dtype='float32'):
    """
    Quickfix to change encoding of the model history for saving in zarr file
    
    :param model_history: model.history from pytorch model (dictionary of pandas DataFrames)
    :param dtype: dtype to convert to
    """
    return {
        key: value.astype(dtype)
        for key, value in model_history.items()
    }


def plot_model_history(title, train, validation, output_path):
    from matplotlib import pyplot as plt
    plt.plot(train, label='train')
    plt.plot(validation, label='validation')
    plt.title(f'Training loss: {title}')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(output_path)
    plt.close()


def check_output(adata, output_type):
    """
    Process data based on output type.
    If more than one output type is given, use the most processed output type: knn > embed > full
    :param adata: integrated anndata object
    :param output_type: string or list of output type
    :return: integrated anndata object with unintegrated anndata in .raw
    """
    if isinstance(output_type, str):
        output_type = [output_type]

    if 'full' in output_type:
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
    elif 'embed' in output_type:
        assert 'X_emb' in adata.obsm
    elif 'knn' in output_type:
        assert 'neighbors' in adata.uns
        assert 'connectivities' in adata.obsp
        assert 'distances' in adata.obsp
    else:
        raise ValueError(f'Invalid output type {output_type}')


def get_hyperparams(hyperparams: dict, model_params: list = None, train_params: list = None):
    """
    Get hyperparameters and training parameters from hyperparameter dictionary
    :param hyperparams: dictionary of hyperparameters
    :param train_params: list of training parameters
    :return: hyperparams, train_params
    """
    if model_params is None:
        model_params = []
        if train_params is not None:
            model_params = [x for x in hyperparams if x not in train_params]
    model_params = {k: v for k, v in hyperparams.items() if k in model_params}
    train_params = {k: v for k, v in hyperparams.items() if k not in model_params}
    return model_params, train_params
