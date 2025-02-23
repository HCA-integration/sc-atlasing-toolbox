import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas.api.types as ptypes
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import tqdm
import traceback
import concurrent.futures
from joblib import Parallel, delayed

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
threads = snakemake.threads


obs_df = read_anndata(input_file, obs='obs').obs
columns = [x for x in cluster_keys + covariates if x in obs_df.columns]
obs_df = obs_df[columns]

def call_plot(covariate, cluster_key, obs_df, output_dir):
    if covariate not in obs_df.columns:
        logging.info(f'Column "{covariate}" not found in obs, skip...')
        return

    try:
        if obs_df[covariate].dtype == 'category':
            score = nmi(obs_df, cluster_key, covariate)
            fig = plot_stacked_bar(
                obs_df,
                category_key=cluster_key,
                covariate_key=covariate,
                count_type='proportion',
                title_suffix=f' (NMI: {score:.3f})',
            )
        elif ptypes.is_numeric_dtype(obs_df[covariate]):
            fig = plot_violin(
                obs_df,
                category_key=cluster_key,
                covariate_key=covariate,
            )
        else:
            raise ValueError(f'Invalid covariate type: {obs_df[covariate].dtype}, must be category or numeric')
        
        fig.savefig(output_dir / f'{cluster_key}--{covariate}.png', bbox_inches='tight')
        
    except Exception as e:
        logging.error(f'Error plotting {cluster_key}--{covariate}: {e}')
        traceback.print_exc()
        return

Parallel(n_jobs=threads, backend='threading')(
    delayed(call_plot)(covariate, cluster_key, obs_df, output_dir)
    for covariate in tqdm(covariates, desc=f'Plotting with {threads} threads')
    for cluster_key in cluster_keys
)

# with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
#     tasks = [
#         (covariate, cluster_key, obs_df, output_dir)
#         for covariate in covariates if covariate in obs_df.columns
#         for cluster_key in cluster_keys
#     ]
#     futures = [executor.submit(call_plot, *task) for task in tasks]
    
#     for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
#         try:
#             future.result()
#         except Exception as e:
#             logging.error(f"Exception occurred: {e}")
#             traceback.print_exc()
