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
file_id = snakemake.wildcards.file_id

logger.info('Read TSV...')
df = pd.read_table(input_file)

if df.shape[0] == 0:
    logger.info('Empty TSV, skip plotting')
    plt.savefig(output_bar)
    plt.savefig(output_violin)
    exit(0)

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
    errorbar='sd',
    dodge=True,
    errwidth=1,
    capsize=.1,
)
g.set(title=f'PCR of covariates for: {dataset} {file_id}')
grouped_by_covariate = df.groupby('covariate', sort=False)
bar_labels = grouped_by_covariate['pcr'].first().round(2).astype(str).str.cat(
    grouped_by_covariate['n_covariates'].first(),
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
plt.grid()
g = sns.violinplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    dodge=False,
)
sns.despine()
g.set(title=f'PCR of covariates for: {dataset} {file_id}')

logger.info('Save violin plot...')
plt.savefig(output_violin, bbox_inches='tight',dpi=300)