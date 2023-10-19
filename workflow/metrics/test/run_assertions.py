import pandas as pd


metrics_columns = ['dataset', 'metric', 'output_type', 'metric_type', 'score']

metrics_df = pd.read_table('test/out/metrics/results/metrics.tsv')

try:
    for column in metrics_columns:
        assert column in metrics_df.columns
except AssertionError:
    raise AssertionError(f'column {column} not found\n{metrics_df}')