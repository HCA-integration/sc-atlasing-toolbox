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

from metrics.utils import anndata_to_mudata, write_metrics
from metrics import metric_map
from utils.io import read_anndata


input_adata = snakemake.input.h5mu
output_file = snakemake.output.metric
dataset = snakemake.wildcards.dataset
file_id = snakemake.wildcards.file_id
metric = snakemake.wildcards.metric
batch_key = snakemake.params.batch_key
label_key = snakemake.params.label_key
# method = snakemake.wildcards.method
# hyperparams = snakemake.wildcards.hyperparams
# lineage_specific = snakemake.wildcards.lineage_specific
# lineage_key = snakemake.wildcards.lineage_key

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']
metric_function = metric_map[metric]

logger.info(f'Read {input_adata} ...')
adata = read_anndata(input_adata)

if 'unintegrated' in snakemake.input.keys():
    input_unintegrated = snakemake.input.unintegrated
    logger.info(f'Read unintegrated data {input_unintegrated}...')
    unintegrated = read_anndata(input_unintegrated)
    # unintegrated = anndata_to_mudata(
    #     unintegrated,
    #     group_key=lineage_key,
    #     # prefix='lineage~'
    # )
else:
    logger.info('Skip unintegrated data...')
    unintegrated = adata

output_types = []
scores = []
# lineages = []
# meta = get_from_adata(adata)
# for output_type in meta['output_types']:
#     if lineage_key == 'None':
#         score = metric_function(
#             mudata[lineage_key],
#             output_type,
#             meta,
#             adata_raw=unintegrated[lineage_key],
#         )
#         scores.append(score)
#         lineages.append('global')
#         output_types.append(output_type)
#     else:
#         for lineage in mudata.mod:
#             score = np.nan
#             if mudata[lineage].obs[meta['label']].nunique() == 1:
#                 continue
#             score = metric_function(
#                 mudata[lineage],
#                 output_type,
#                 meta,
#                 adata_raw=unintegrated[lineage],
#             )
#             scores.append(score)
#             lineages.append(lineage)
#             output_types.append(output_type)

for output_type in adata.uns['output_types']:
    score = metric_function(
        adata,
        output_type,
        batch_key=batch_key,
        label_key=label_key,
        adata_raw=unintegrated,
    )
    scores.append(score)
    output_types.append(output_type)


write_metrics(
    scores=scores,
    output_types=output_types,
    metric=metric,
    metric_type=metric_type,
    batch=batch_key,
    label=label_key,
    file_id=file_id,
    dataset=dataset,
    filename=output_file
)
