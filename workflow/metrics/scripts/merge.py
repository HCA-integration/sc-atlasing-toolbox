import re
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
ex_columns = set()

# for col in expanded_file_ids.columns:
#     if expanded_file_ids[col].str.contains('=').any():
#         exp_col = expanded_file_ids[col].str.split('=', expand=True)
#     else:
#         exp_col = expanded_file_ids[col].str.rsplit(':', n=1, expand=True)
#     non_empty = exp_col[0].apply(lambda x: isinstance(x, str))
#     col_name = exp_col[non_empty].iloc[0,0]
#     ex_columns.append(col_name)
#     # reassign value to column
#     expanded_file_ids[col] = exp_col.iloc[:,-1]
# # rename expanded columns
# expanded_file_ids.columns = ex_columns
# # sort columns by name
# ex_columns.sort()
# expanded_file_ids = expanded_file_ids[ex_columns]
# # merge extra columns to metrics
# if 'file_id' in expanded_file_ids.columns:
#     del metrics_df['file_id']
# metrics_df = pd.concat([metrics_df, expanded_file_ids], axis=1)

for row, _dict in expanded_file_ids.to_dict('index').items():
    for _, value in _dict.items():
        if value is None:
            continue
        splits = re.split(r'[:=]', value, maxsplit=1)
        if len(splits) == 1:
            continue
        key, value = splits
        metrics_df.loc[row, key] = value
        ex_columns.add(key)

# save files
with open(extra_columns, 'w') as f:
    for col in sorted(ex_columns):
        f.write(f'{col}\n')
metrics_df.to_csv(out_tsv, sep='\t', index=False)

print(metrics_df[['metric', 'output_type', 'metric_type', 'score', 's', 'h:m:s']])
