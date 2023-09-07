# Preprocessing

This module will further process datasets following the standard preprocessing
steps: normalize counts, calculate highly variable genes, PCA, neighbors, and
the additional dimentionality reduction (UMAP).

## General input

It must include:

* `.X` with raw counts or specify in config file as `raw_counts`.
* `.obs` batch and lineage columns if user wishes to run
`highly_variable_genes` with them.

Also, it can take the anndata.Anndata file the `load_data` workflow outputs.

## Steps' explanation

### Normalize

Transforms the data using normalize_total and then log1p. It also makes sure it is a sparse matrix.

**Output**

- `.X` normalised and log1 transformed counts.

**Parameters**

- raw_counts: H5AD file to raw counts.

### Highly variable genes selection

Highly variable genes are calculated after using `sc.pp.filter_genes(min_cells=1)` and the `.var` in the unfiltered object is updated.

**Output**

- `.var["highly_variable"]` boolean column with highly variable gene status.

**Parameters**

- batch: batch column to account when getting the HVGs.
- lineage: if you want to calculate lineage specific genes alone
  or combined with batch.
- args: more to pass to the `scanpy.pp.highly_variable_genes` function.

### PCA

If `scale=True` scaling of `.X` will be performed. Then the PCA is
calculated using the highly varaible genes.

**Output**

- `.obsm["X_pca"]` PCA embedding

**Parameters**

- scale: wether or not to scale counts before calculating
the PCA embedding.

### Neighbors

It will attempt to use the RAPIDS[^1] implementation but will default, if it fails, to the UMAP implementation
[arXiv:1802.03426v3](https://arxiv.org/abs/1802.03426v3).

**Output**

- `.obsp["neighbors"]["distances"]` neighborhood graph
- `.obsp["neighbors"]["connectivities"]` neighborhood graph

**Parameters**

- args: extra arguments for the `scanpy.pp.neighbors` function.

### UMAP

UMAP dimensionality reduction is calculated from the PCA output. It can also use the RAPIDS[^1] implementation.

**Output**

- `.obsm["X_umap"]` UMAP embedding; it becoms .obsm["X_umap_{key}"] if multiple `neighbors_key` are given in the config file.

## General output

anndata.Anndata with:

- `.uns["preprocessing"]` adds 'highly_variable_genes', 'normalization',
'scaled' and 'log-transformed' info.

Also, the components generated in each step as described above can be
combined  by instructing in the config file like so:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
    assemble:
      - highly_variable_genes
      - pca
      - neighbors
```

---

[^1] On RAPIDS implementation: whether it is used or not depends on 'os'
in the config file.
