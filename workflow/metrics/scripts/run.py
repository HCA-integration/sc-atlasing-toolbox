import pandas as pd
import mudata as mu
from metrics.utils import write_metrics, get_from_adata
from metrics import metric_map

input_adata = snakemake.input.h5mu
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

print(f'read {input_adata} ...')
mudata = mu.read(input_adata)
meta = get_from_adata(mudata)

output_types = []
scores = []
lineages = []
for output_type in meta['output_types']:
    if lineage_key == 'None':
        score = metric_function(
            mudata[lineage_key],
            output_type,
            meta
        )
        scores.append(score)
        lineages.append('global')
        output_types.append(output_type)
    else:
        for lineage in mudata.mod:
            score = metric_function(
                mudata[lineage],
                output_type,
                meta
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
