import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import rcParams
# figure size in inches
rcParams['figure.figsize'] = 12, 9


input_file = snakemake.input.tsv
output_bar = snakemake.output.barplot
output_violin = snakemake.output.violinplot
dataset = snakemake.wildcards.dataset

logger.info('Read TSV...')
df = pd.read_table(input_file)
df = df.sort_values(
    ['pcr', 'n_covariates', 'covariate'],
    ascending=[False, True, False]
)

logger.info('Barplot...')
g = sns.barplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    dodge=True,
    errwidth=0,
    # errorbar='sd',
)
g.set(title=f'PCR of covariates for: {dataset}')
bar_labels = df.groupby('covariate')['pcr'].mean().round(2).astype(str).str.cat(
    df.groupby('covariate')['n_covariates'].first(),
    sep=', '
)
g.bar_label(
    g.containers[0],
    labels=bar_labels,
    padding=10
)
plt.xticks(rotation=90)
sns.despine()

logger.info('Save barplot...')
plt.savefig(output_bar, bbox_inches='tight',dpi=300)

logging.info('Violin plot...')
plt.clf()
g = sns.violinplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    dodge=False,
)
plt.grid()

logger.info('Save barplot...')
plt.savefig(output_violin, bbox_inches='tight',dpi=300)