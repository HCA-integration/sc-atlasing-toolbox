import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import pandas as pd

from utils.io import read_anndata, write_zarr_linked
from qc_utils import get_thresholds, apply_thresholds


def thresholds_to_df(df, wildcards, qc_type=None, **kwargs):
    if qc_type is None: 
        if 'user_thresholds' in kwargs:
            qc_type = 'user'
            if 'autoqc_thresholds' in kwargs:
                qc_type = 'updated'
        elif 'autoqc_thresholds' in kwargs:
            qc_type = 'sctk_autoqc'
        else:
            qc_type = str(qc_type)
    
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
output_tsv = snakemake.output.tsv
output_qc_stats = snakemake.output.qc_stats

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', uns='uns')

threshold_keys = ['n_counts', 'n_genes', 'percent_mito']
user_thresholds = snakemake.params.get('thresholds')
alternative_thresholds = snakemake.params.get('alternative_thresholds')
autoqc_thresholds = adata.uns.get('scautoqc_ranges')

# Calculate threshold stats
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
        # alternative thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            qc_type='alternative',
            threshold_keys=threshold_keys,
            user_thresholds=alternative_thresholds,
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
adata.uns['qc'] = df

logging.info(f'Write thresholds to {output_tsv}...')
df.to_csv(output_tsv, sep='\t', index=False)

# add QC status column to .obs with 'passed', 'failed', and 'ambiguous'
apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        user_thresholds=user_thresholds,
        autoqc_thresholds=autoqc_thresholds,
    ),
    threshold_keys=threshold_keys,
    column_name='user_qc_status',
)
apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        user_thresholds=alternative_thresholds,
        autoqc_thresholds=autoqc_thresholds,
    ),
    threshold_keys=threshold_keys,
    column_name='alternative_qc_status',
)

# get 'passed' if user_qc_status and alternative_qc_status are both True
user_status = adata.obs['user_qc_status']
alt_status = adata.obs['alternative_qc_status']
adata.obs['qc_status'] = 'ambiguous'
adata.obs.loc[user_status & alt_status, 'qc_status'] = 'passed'
adata.obs.loc[~(user_status | alt_status), 'qc_status'] = 'failed'

qc_status_counts = pd.Categorical(
    adata.obs['qc_status'],
    categories=['passed', 'failed', 'ambiguous']
).value_counts()
qc_status_counts = pd.DataFrame(qc_status_counts).T
qc_status_counts['file_id'] = snakemake.wildcards.file_id
qc_status_counts = qc_status_counts[['file_id', 'passed', 'failed', 'ambiguous']]

logging.info(f'Write QC stats to {output_qc_stats}...')
qc_status_counts.to_csv(output_qc_stats, sep='\t', index=False)

logging.info(f'Write {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obs', 'uns']
)