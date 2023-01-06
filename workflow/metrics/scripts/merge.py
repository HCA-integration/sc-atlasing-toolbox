import pandas as pd
import warnings

warnings.filterwarnings('ignore')

input_metrics = snakemake.input.metrics
input_benchmark = snakemake.input.benchmark
out_tsv = snakemake.output.tsv
wildcards = snakemake.wildcards
expanded_wildcards = snakemake.params['wildcards']
expanded_wildcards = pd.DataFrame(expanded_wildcards)
facet_var = expanded_wildcards.columns[0]

metrics_df = pd.concat(
    [pd.read_table(file) for file in input_metrics],
    ignore_index=True,
).reset_index(drop=True)

print(metrics_df)

benchmark_df = pd.concat(
    [pd.read_table(file) for file in input_benchmark],
    ignore_index=True,
).reset_index(drop=True)

print(benchmark_df)

benchmark_df = pd.concat([expanded_wildcards, benchmark_df], axis=1)

metrics_df = metrics_df.merge(benchmark_df)
metrics_df.to_csv(out_tsv, sep='\t', index=False)

print(metrics_df[['metric', 'method', 'output_type', 'metric_type', 'score', 's', 'h:m:s']])
