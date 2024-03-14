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


input_file = snakemake.input[0]
output_file = snakemake.output[0]
dataset = snakemake.wildcards.dataset
file_id = snakemake.wildcards.file_id
metric = snakemake.wildcards.metric
params = snakemake.params
batch_key = params.batch_key
label_key = params.label_key

# metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = params.get('metric_type')
assert metric_type in ['batch_correction', 'bio_conservation'], f'Unknown metric_type: {metric_type}'
metric_output_types = params.get('output_types', [])
comparison = params.get('comparison', False)
metric_function = metric_map.get(metric)

uns = read_anndata(input_file, uns='uns').uns
data_output_type = uns.get('output_type', 'full') # Same as in prepare.py
if data_output_type not in metric_output_types:
    logging.info(
        f'Skip metric={metric} for data output type={data_output_type}\nmetrics output types={metric_output_types}'
    )
    write_metrics(
        scores=[np.nan],
        output_types=[data_output_type],
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
    obsp='obsp',
    var='var',
    uns='uns',
)
if 'embed' in metric_output_types:
    kwargs |= {'obsm': 'obsm'}
if 'full' in metric_output_types:
    kwargs |= {'X': 'X'}
if comparison:
    kwargs |= {'raw': 'raw', 'varm': 'varm', 'obsm': 'obsm'}

adata = read_anndata(input_file, **kwargs)
if comparison:
    adata_raw = adata.raw.to_adata()
    adata_raw.obs = adata.obs
    adata_raw.var = adata.var
    adata_raw.uns = adata.uns
    adata_raw.varm = adata.varm
else:
    adata_raw = None

data_output_type = adata.uns.get('output_type', 'full')
logger.info(f'Run metric {metric} for {data_output_type}...')
adata.obs[label_key] = adata.obs[label_key].astype(str).fillna('NA').astype('category')
score = metric_function(
    adata,
    data_output_type,
    batch_key=batch_key,
    label_key=label_key,
    adata_raw=adata_raw,
)

write_metrics(
    scores=[score],
    output_types=[data_output_type],
    metric=metric,
    metric_type=metric_type,
    batch=batch_key,
    label=label_key,
    file_id=file_id,
    dataset=dataset,
    filename=output_file,
    **adata.uns.get('wildcards', {}),
)
