import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd


def check_same_categories(df, columns):
    first_categories = df[columns[0]].cat.categories
    if not all(df[col].cat.categories.equals(first_categories) for col in columns):
        logging.warn(f"Columns {columns} have different categories, majority voting might not work correctly")


def get_majority_reference(df, reference_key, query_key, **kwargs):
    map_majority = pd.crosstab(df[reference_key], df[query_key], **kwargs).idxmax(axis=0)
    return pd.Categorical(
        df[query_key].map(map_majority),
        categories=map_majority.dropna().unique(),
    )

def get_majority_consensus(df, columns, new_key='majority_consensus'):
    from scipy import stats
    from functools import reduce
    
    agreement_col = f'{new_key}_agreement'
    low_agreement_col = f'{new_key}_low_agreement'
    n_vote_columns = len(columns)
    
    check_same_categories(df, columns)
    
    categories = reduce(pd.Index.union, [df[col].cat.categories.dropna() for col in columns])
    cat_dtype = pd.CategoricalDtype(categories=categories)
    df_cat = df[columns].astype(cat_dtype).apply(lambda x: x.cat.codes)
    
    logging.info(f'Compute majority vote across {n_vote_columns} columns...')
    majority_votes = stats.mode(df_cat, axis=1, keepdims=False)
    
    logging.info('Assign majority votes to new column...')
    df[new_key] = pd.Categorical.from_codes(majority_votes.mode, categories=categories)
    df[agreement_col] = majority_votes.count / n_vote_columns
    
    min_agreement = 1 / n_vote_columns
    df.loc[df[agreement_col] <= min_agreement, new_key] = float('nan')
    # set to 0, because no agreement when all assignments different
    df.loc[df[agreement_col] <= min_agreement, agreement_col] = 0
    df[low_agreement_col] = df[agreement_col] <= 0.5
    
    counts = df.value_counts(
        subset=[new_key, agreement_col],
        sort=False,
        dropna=False
    ).reset_index(name='count')
    counts['ratio'] = counts['count'] / counts.groupby(new_key, observed=True)['count'].transform('sum')
    print(counts, flush=True)
    counts.to_csv(output_plots / 'majority_consensus_counts.tsv', sep='\t', index=False)
    
    return df[[new_key, agreement_col, low_agreement_col]]
    

if __name__ == '__main__':
    import seaborn as sns
    from matplotlib import pyplot as plt
    
    from utils.io import read_anndata, write_zarr_linked

    input_file = snakemake.input.zarr
    output_file = snakemake.output.zarr
    output_plots = Path(snakemake.output.plots)
    output_plots.mkdir(parents=True, exist_ok=True)
    majority_reference = snakemake.params.get('majority_reference')
    majority_consensus = snakemake.params.get('majority_consensus')

    logging.info('Read adata...')
    adata = read_anndata(input_file, obs='obs')

    if majority_reference is not None:
        adata.obs['majority_reference'] = get_majority_reference(
            adata.obs,
            reference_key=majority_reference['reference_key'],
            query_key=majority_reference['query_key'],
            **majority_reference.get('crosstab_kwargs', {})
        )

    if majority_consensus is not None:
        new_key = 'majority_consensus'
        maj_df = get_majority_consensus(
            adata.obs,
            columns=majority_consensus['columns'],
            new_key=new_key,
            **majority_consensus.get('mode_kwargs', {})
        )
        adata.obs[maj_df.columns] = maj_df
        
        plot_df = pd.crosstab(
            maj_df[new_key],
            maj_df[f'{new_key}_low_agreement'],
            dropna=False,
            normalize='index'
        )
        
        if True not in plot_df.columns:
            plot_df[True] = 0
        elif False not in plot_df.columns:
            plot_df[False] = 0
        
        ax = plot_df.sort_values(True).plot(
            kind='barh',
            stacked=True,
            color=sns.color_palette("colorblind").as_hex(),
            figsize=(5, 7),
        )
        ax.legend(title='Low agreement', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(False)    
        plt.savefig(output_plots / 'majority_consensus_frac.png', bbox_inches='tight')
                
    logging.info(f'Write to {output_file}...')
    write_zarr_linked(
        adata=adata,
        in_dir=input_file,
        out_dir=output_file,
        files_to_keep=['obs']
    )