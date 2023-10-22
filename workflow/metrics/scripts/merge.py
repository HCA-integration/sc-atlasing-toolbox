import pandas as pd
import warnings

warnings.filterwarnings('ignore')

input_metrics = snakemake.input.metrics
input_benchmark = snakemake.input.benchmark
out_tsv = snakemake.output.tsv
extra_columns = snakemake.output.extra_columns
expanded_wildcards = snakemake.params['wildcards']
# expanded_wildcards = pd.DataFrame(expanded_wildcards)

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

metrics_df = metrics_df.merge(benchmark_df).drop_duplicates()

# parse any extra entries in file_id
expanded_file_ids = metrics_df['file_id'].str.split('--', expand=True)
ex_columns = []
prefix = ''

for col in expanded_file_ids.columns:
    exp_col = expanded_file_ids[col].str.split('=', expand=True)
    expanded_file_ids[col] = exp_col.iloc[:,-1]
    col_name = exp_col.iloc[0,0]
    if ':' in col_name:
        prefix = col_name.split(':')[0]
        col_name = col_name.replace(':', '_')
    else:
        col_name = f'{prefix}_{col_name}'
    ex_columns.append(col_name)

if len(ex_columns) == 1:
    ex_columns = ['file_id']

expanded_file_ids.columns = ex_columns
metrics_df = pd.concat([metrics_df, expanded_file_ids], axis=1)

# save files
with open(extra_columns, 'w') as f:
    for col in ex_columns:
        f.write(f'{col}\n')
metrics_df.to_csv(out_tsv, sep='\t', index=False)

print(metrics_df[['metric', 'output_type', 'metric_type', 'score', 's', 'h:m:s']])
