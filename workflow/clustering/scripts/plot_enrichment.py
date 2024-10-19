import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas.api.types as ptypes
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import tqdm

from utils.io import read_anndata
from utils.plots import plot_stacked_bar, plot_violin


def nmi(df, x, y):
    from sklearn.metrics.cluster import normalized_mutual_info_score

    df = df[[x, y]]
    df = df[df[x].notna() & df[y].notna()]
    return normalized_mutual_info_score(df[x], df[y])


input_file = snakemake.input.zarr
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(exist_ok=True)
cluster_keys = snakemake.params.get('cluster_keys')
covariates = snakemake.params.get('covariates')


obs_df = read_anndata(input_file, obs='obs').obs

for covariate in tqdm(covariates):
    if covariate not in obs_df.columns:
        logging.info(f'Column "{covariate}" not found in obs, skip...')
        continue
    for cluster_key in cluster_keys:
        try:
            if obs_df[covariate].dtype == 'category':
                score = nmi(obs_df, cluster_key, covariate)
                plot_stacked_bar(
                    obs_df,
                    category_key=cluster_key,
                    covariate_key=covariate,
                    count_type='proportion',
                    title_suffix=f' (NMI: {score:.3f})',
                )
            elif ptypes.is_numeric_dtype(obs_df[covariate]):
                # sns.violinplot(x=covariate, y=cluster_key, data=obs_df)
                plot_violin(
                    obs_df,
                    category_key=cluster_key,
                    covariate_key=covariate,
                )
            else:
                raise ValueError(f'Invalid covariate type: {obs_df[covariate].dtype}')
        except TypeError as e:
            logging.error(f'Error plotting {cluster_key}--{covariate}: {e}')
            continue
        plt.savefig(output_dir / f'{cluster_key}--{covariate}.png', bbox_inches='tight')
        plt.show()