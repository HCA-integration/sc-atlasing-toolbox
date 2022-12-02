import pandas as pd
import glob

# celltypist
celltypist = glob.glob('test/out/label_transfer/celltypist/*/*.tsv')
if len(celltypist) == 0:
    print('No files found for celltypist')

for file in celltypist:
    print(f'read {file}...')
    predictions = pd.read_table(file)
    celltypist_columns = [x for x in predictions.columns if x.startswith('celltypist')]

    try:
        assert len(celltypist_columns) == 4
    except AssertionError:
        print('Error: Incorrect number of columns')
        print(predictions)
