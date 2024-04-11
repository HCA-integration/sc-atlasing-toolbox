import logging
import pandas as pd
import cellhint

from utils.io import read_anndata

logger = logging.getLogger(__name__)

input_file = snakemake.input[0]
output_file = snakemake.output.reannotation
output_reannotation = snakemake.output.reannotation

print('read...')
df = pd.read_table(input_file)
mapping = {k: str(i) for i, k in enumerate(df['reannotation'].unique())}
df['reannotation_index'] = df['reannotation'].map(mapping)

print(df)

df.to_csv(output_reannotation, sep='\t', index=False)
