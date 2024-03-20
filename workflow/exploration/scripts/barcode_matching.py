import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from anndata.experimental import read_elem
import zarr


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


input_file = snakemake.input.zarr
output_png = snakemake.output.png

z = zarr.open(input_file)
obs = read_elem(z["obs"])

# clean obsnames
obs_names = pd.Series(obs.index).str.split('-', n=1, expand=True)[0]

# Select random 100 samples to plot
np.random.seed(42)
size = np.min((100, obs['sample'].nunique()))
samples = np.random.choice(obs['sample'].unique(), size=size, replace=False).tolist()
barcodes = [obs.query(f'sample == "{sample}"').index.tolist() for sample in samples]

UpSetFromLists(
    listOflist=barcodes,
    labels=samples,
    verbose=False,
).plot()
plt.savefig(output_png)
