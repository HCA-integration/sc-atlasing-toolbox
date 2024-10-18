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


def nmi(df, x, y):
    from sklearn.metrics.cluster import normalized_mutual_info_score

    df = df[[x, y]]
    df = df[df[x].notna() & df[y].notna()]
    return normalized_mutual_info_score(df[x], df[y])


def plot_stacked_bar(df, x, y, count_type='count', palette='tab20', title=''):
    colors = palette
    if palette is not None:
        colors = plt.get_cmap(palette).colors
    df = df[[x, y]].copy()
    df[x] = df[x].astype(str)
    counts_df = df[[x, y]].value_counts().to_frame()
    if count_type == 'proportion':
        counts_df = counts_df / counts_df.groupby(x, observed=True).transform('sum')
        counts_df.columns = ['proportion']
    counts_df = counts_df.reset_index(drop=False)
    
    num_categories = counts_df[x].nunique()
    plt.figure(figsize=(8, max(5, num_categories * 0.3)))
    ax = counts_df.pivot(
        index=x,
        columns=y,
        values=count_type
    ).fillna(0) \
    .plot.barh(
        stacked=True,
        color=colors,
        ax=plt.gca()
    )
        
    if df[y].nunique() < 30:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        ax.legend_ = None
    title = f'{count_type.capitalize()} of {y} by {x}{title}'
    ax.set_title(title)


def plot_violin(df, cluster_key, covariate_key):
    clusters = sorted(df[cluster_key].unique())
    grouped_data = [
        df[df[cluster_key] == cluster][covariate_key]
        for cluster in clusters
    ]
    num_categories = len(grouped_data)
    plt.figure(figsize=(8, max(5, num_categories * 0.3)))
    plt.violinplot(
        grouped_data,
        showmeans=False,
        showmedians=True,
        vert=False,
    )
    plt.yticks(ticks=range(1, len(grouped_data) + 1), labels=clusters)
    plt.title(f'Distribution of {covariate_key} by {cluster_key}')
    plt.xlabel(covariate_key)
    plt.ylabel(cluster_key)


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
                    cluster_key,
                    covariate,
                    count_type='proportion',
                    title=f' (NMI: {score:.3f})',
                )
            elif ptypes.is_numeric_dtype(obs_df[covariate]):
                # sns.violinplot(x=covariate, y=cluster_key, data=obs_df)
                plot_violin(obs_df, cluster_key, covariate)
            else:
                raise ValueError(f'Invalid covariate type: {obs_df[covariate].dtype}')
        except TypeError as e:
            logging.error(f'Error plotting {cluster_key}--{covariate}: {e}')
            continue
        plt.savefig(output_dir / f'{cluster_key}--{covariate}.png', bbox_inches='tight')
        plt.show()