import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt


def UpSetFromLists(listOflist, labels, max_overlap=5, min_overlap_size=50, verbose=False, **kwargs):
    """
    code taken and adapted from https://github.com/brianpenghe/python-genomics/blob/master/pandasPlus.py
    :param listOflist:
    :param labels:
    :param kwargs:
    :return:
    """
    import numpy as np
    from upsetplot import UpSet

    listall = list(set([j for i in listOflist for j in i]))
    temp = pd.Series(listall, index=listall)

    temp2 = pd.concat([temp.isin(i) for i in listOflist + [temp]], axis=1)
    temp2.columns = labels + ['all']
    upset_counts = temp2.value_counts()

    upset_counts = upset_counts[[sum(tup) <= max_overlap for tup in upset_counts.index]]
    min_overlap_size = np.min([upset_counts.max(), min_overlap_size])
    upset_counts = upset_counts[upset_counts >= min_overlap_size]

    if verbose:
        print('index', upset_counts.index)
        print('# duplicated', upset_counts.index.duplicated().sum())
        print('value_counts', upset_counts.reset_index(drop=True))

    return UpSet(upset_counts, **kwargs)


input_h5ad = snakemake.input.h5ad
output_tsv = snakemake.output.tsv

adata = sc.read(input_h5ad)

# clean obsnames
adata.obs_names = pd.Series(adata.obs_names).str.split('-', n=1, expand=True)[0]
adata.obs['sample'].reset_index().to_csv(output_tsv, sep='\t', index=False)

# samples = adata.obs['sample'].unique().to_list()
# barcodes = [adata[adata.obs['sample'] == sample].obs_names.tolist() for sample in samples]
#
# UpSetFromLists(
#     listOflist=barcodes,
#     labels=samples,
#     verbose=True
# ).plot()
# plt.savefig(output_png)
