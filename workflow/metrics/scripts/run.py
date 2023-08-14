import numpy as np
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO)

from metrics.utils import anndata_to_mudata, write_metrics, get_from_adata
from metrics import metric_map
from utils.io import read_anndata_or_mudata


input_adata = snakemake.input.h5mu
input_unintegrated = snakemake.input.unintegrated
output_file = snakemake.output.metric
metric = snakemake.wildcards.metric
method = snakemake.wildcards.method
dataset = snakemake.wildcards.dataset
hyperparams = snakemake.wildcards.hyperparams
lineage_specific = snakemake.wildcards.lineage_specific
lineage_key = snakemake.wildcards.lineage_key

metrics_meta = pd.read_table(snakemake.input.metrics_meta, index_col='metric')
metric_type = metrics_meta.loc[metric]['metric_type']
metric_function = metric_map[metric]

logging.info(f'Read {input_adata} ...')
mudata = read_anndata_or_mudata(input_adata)
meta = get_from_adata(mudata)

if metrics_meta.query(f'metric == "{metric}"')['comparison'].all():
    logging.info(f'Read unintegrated data {input_unintegrated}...')
    unintegrated = read_anndata_or_mudata(input_unintegrated)
    unintegrated = anndata_to_mudata(
        unintegrated,
        group_key=lineage_key,
        # prefix='lineage~'
    )
else:
    logging.info('Skip unintegrated data...')
    unintegrated = mudata

output_types = []
scores = []
lineages = []
for output_type in meta['output_types']:
    if lineage_key == 'None':
        score = metric_function(
            mudata[lineage_key],
            output_type,
            meta,
            adata_raw=unintegrated[lineage_key],
        )
        scores.append(score)
        lineages.append('global')
        output_types.append(output_type)
    else:
        for lineage in mudata.mod:
            score = np.nan
            if mudata[lineage].obs[meta['label']].nunique() == 1:
                continue
            score = metric_function(
                mudata[lineage],
                output_type,
                meta,
                adata_raw=unintegrated[lineage],
            )
            scores.append(score)
            lineages.append(lineage)
            output_types.append(output_type)

write_metrics(
    scores=scores,
    output_types=output_types,
    lineages=lineages,
    lineage_specific=lineage_specific,
    lineage_key=lineage_key,  # TODO: better name for per lineage specific evaluation
    metric=metric,
    metric_type=metric_type,
    method=method,
    hyperparams=hyperparams,
    dataset=dataset,
    filename=output_file
)
