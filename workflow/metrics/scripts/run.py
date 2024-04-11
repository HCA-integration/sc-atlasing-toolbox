import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('Run metric')
logger.setLevel(logging.INFO)

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    logger.info('no hardware acceleration for sklearn')

from metrics import metric_map
from metrics.utils import write_metrics
from utils.io import read_anndata
from utils.accessors import subset_hvg


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

dataset = wildcards.dataset
file_id = wildcards.file_id
batch_key = wildcards.batch
label_key = wildcards.label
metric = wildcards.metric

metric_type = params.get('metric_type')
assert metric_type in ['batch_correction', 'bio_conservation'], f'Unknown metric_type: {metric_type}'
allowed_output_types = params.get('output_types')
input_type = params.get('input_type')
comparison = params.get('comparison', False)
metric_function = metric_map.get(metric, ValueError(f'No function for metric: {metric}'))

uns = read_anndata(input_file, uns='uns').uns
output_type = uns.get('output_type', 'full') # Same as in prepare.py

if output_type not in allowed_output_types:
    logging.info(
        f'Skip metric={metric} for data output type={output_type}\n'
        f'allowed output types={allowed_output_types}'
    )
    write_metrics(
        scores=[np.nan],
        output_types=[output_type],
        metric=metric,
        metric_type=metric_type,
        batch=batch_key,
        label=label_key,
        file_id=file_id,
        dataset=dataset,
        filename=output_file,
        **uns.get('wildcards', {}),
        )
    exit(0)

logger.info(f'Read {input_file} ...')
kwargs = dict(
    obs='obs',
    uns='uns',
)
if input_type == 'knn':
    kwargs |= {'obsp': 'obsp'}
if input_type == 'embed':
    kwargs |= {'obsm': 'obsm'}
if input_type == 'full':
    kwargs |= {'X': 'X', 'var': 'var'}

adata = read_anndata(input_file, **kwargs)
print(adata, flush=True)
adata_raw = None

if comparison:
    kwargs = {'obs': 'obs', 'var': 'raw/var', 'varm': 'raw/varm', 'X': 'raw/X'}
    adata_raw = read_anndata(input_file, **kwargs, dask=True, backed=True)
    print('adata_raw')
    print(adata_raw, flush=True)
    
    adata_raw, _ = subset_hvg(adata_raw, var_column='highly_variable', to_memory=[], compute_dask=False)
    adata_raw.obs = adata.obs
    adata_raw.obsm = adata.obsm
    adata_raw.uns = adata.uns
    adata_raw.varm = adata.varm
else:
    adata_raw = None

logger.info(f'Run metric {metric} for {output_type}...')
adata.obs[label_key] = adata.obs[label_key].astype(str).fillna('NA').astype('category')
score = metric_function(
    adata,
    output_type,
    batch_key=batch_key,
    label_key=label_key,
    adata_raw=adata_raw,
)

write_metrics(
    scores=[score],
    output_types=[output_type],
    metric=metric,
    metric_type=metric_type,
    batch=batch_key,
    label=label_key,
    file_id=file_id,
    dataset=dataset,
    filename=output_file,
    **adata.uns.get('wildcards', {}),
)
