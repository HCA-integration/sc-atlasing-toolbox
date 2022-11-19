import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from upsetplot import UpSet


def UpSetFromLists(listOflist, labels, **kwargs):
    """
    code taken and adapted from https://github.com/brianpenghe/python-genomics/blob/master/pandasPlus.py
    :param listOflist:
    :param labels:
    :param kwargs:
    :return:
    """
    listall = list(set([j for i in listOflist for j in i]))
    temp = pd.Series(listall, index=listall)
    temp2 = pd.concat([temp.isin(i) for i in listOflist + [temp]], axis=1)
    temp2.columns = labels + ['all']
    temp2 = temp2.set_index(labels)
    return UpSet(temp2, subset_size='count', intersection_plot_elements=3, **kwargs)


input_h5ad = snakemake.input.h5ad
output_png = snakemake.output.png

adata = sc.read(input_h5ad)
# clean obsnames
adata.obs_names = pd.Series(adata.obs_names).str.split('-', n=1, expand=True)[0]

samples = adata.obs['sample'].unique().to_list()
barcodes = [adata[adata.obs['sample'] == sample].obs_names.tolist() for sample in samples]

UpSetFromLists(
    listOflist=barcodes,
    labels=samples,
).plot()
plt.savefig(output_png)
