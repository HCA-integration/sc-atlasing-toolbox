def by_sample(adata, n_cell_max, sample_key, random_state=42):
    """
    Randomly subset complete samples until the maximum number of cells is reached
    """
    samples = []
    n_cells = 0

    shuffled_samples = adata.obs[sample_key].value_counts().sample(frac=1, random_state=random_state)
    
    for sample, count in shuffled_samples.items():
        if len(samples) > 0 and n_cells > n_cell_max:
            break
        n_cells += count
        samples.append(sample)

    return adata.obs[sample_key].isin(samples)


def within_sample(adata, n_cell_max, sample_key, random_state=42):
    """
    Subset to n_cell_max / n_samples random cells per sample
    """
    n_samples = adata.obs[sample_key].nunique()
    n_cells_per_sample = int(n_cell_max / n_samples)
    adata.obs['subset'] = False

    for sample in adata.obs[sample_key].unique():
        # subset to sample
        sample_mask = adata.obs[sample_key] == sample
        # turn mask into subset index
        sample_mask = sample_mask[sample_mask].index.to_series()
        
        # subset to n_cells_per_sample random cells
        n_cells = min(len(sample_mask), n_cells_per_sample)
        subset_sample_mask = sample_mask.sample(n_cells, random_state=random_state)

        # flag cells to be included in subset        
        adata.obs.loc[subset_sample_mask, 'subset'] = True

    return adata.obs['subset']


def scarf_TopACeDo(adata, n_cell_max, sample_key):
    pass


SUBSET_MAP = {
    'by_sample': by_sample,
    'within_sample': within_sample,
    'scarf_TopACeDo': scarf_TopACeDo,
}