import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd

from utils.io import read_anndata
from qc_utils import get_thresholds


def thresholds_to_df(df, wildcards, **kwargs):
    qc_type = 'None'
    if 'user_thresholds' in kwargs:
        kwargs['user_key'] = wildcards.file_id
        qc_type = 'user'
        if 'autoqc_thresholds' in kwargs:
            qc_type = 'updated'
    elif 'autoqc_thresholds' in kwargs:
        qc_type = 'sctk_autoqc'
    
    # get thresholds
    thresholds = get_thresholds(**kwargs, transform=False)
    
    # sort threshold columns
    threshold_keys = sorted(list(thresholds.keys()), reverse=True)
    thresholds = {key: thresholds[key] for key in threshold_keys}
    
    if df.shape[0] == 0:
        thresholds = dict(
            dataset=wildcards.dataset,
            file_id=wildcards.file_id,
            threshold_type=qc_type,
            passed_frac=0,
            removed_frac=0,
            n_passed=0,
            n_removed=0,
            n_total=0,
        ) | thresholds
        return pd.DataFrame(thresholds, index=[0])
    
    # apply thresholds
    df['passed_qc'] = True
    for key in kwargs['threshold_keys']:
        df['passed_qc'] = df['passed_qc'] & df[key].between(thresholds[f'{key}_min'], thresholds[f'{key}_max'])
    passed_qc = pd.Categorical(df['passed_qc'], categories=[True, False]).value_counts()
    
    thresholds = dict(
        dataset=wildcards.dataset,
        file_id=wildcards.file_id,
        threshold_type=qc_type,
        passed_frac=passed_qc[True] / passed_qc.sum(),
        removed_frac=passed_qc[False] / passed_qc.sum(),
        n_passed=passed_qc[True],
        n_removed=passed_qc[False],
        n_total=passed_qc.sum(),
    ) | thresholds
    
    return pd.DataFrame(thresholds, index=[0])


input_file = snakemake.input[0]
output_file = snakemake.output[0]

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', uns='uns')

threshold_keys = ['n_counts', 'n_genes', 'percent_mito']
user_thresholds = snakemake.params.get('thresholds')
autoqc_thresholds = adata.uns.get('scautoqc_ranges')

df = pd.concat(
    [
        # autoqc thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            autoqc_thresholds=autoqc_thresholds,
        ),
        # user thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            user_thresholds=user_thresholds,
        ), 
        # updated thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            autoqc_thresholds=autoqc_thresholds,
            user_thresholds=user_thresholds,
        ), 
    ]
)
print(df, flush=True)
df.to_csv(output_file, sep='\t', index=False)
