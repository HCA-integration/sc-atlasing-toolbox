import pandas as pd


metrics_columns = ['metric', 'method', 'output_type', 'metric_type', 'score']

metrics_df = pd.read_table('test/out/metrics/metrics.tsv')

for column in metrics_df.columns:
    assert column in metrics_columns

for column in metrics_columns:
    assert column in metrics_df.columns