# Preprocessing

This module will further process datasets following the standard preprocessing
steps: normalize counts, calculate highly variable genes, PCA, neighbors, and
the additional dimentionality reduction (UMAP).

## Input

It will take the output from the anndata.Anndata file the `load_data` workflow
outputs. It must include:

* `.X` with raw counts
* `.obs` containing all the columns that are defined by `load_data` workflow
* `.var` containing all the columns that are defined by `load_data` workflow

## Output

anndata.Anndata with:

* `.X` normalised and log1 transformed counts from rule `normalize`.
* `.var["highly_variable"]` highly variable gene status from rule
`highly_variable_genes`.
* `.obsm["X_pca"]` PCA embedding from rule `pca`.
* `.obsp["neighbors"]["distances"]` neighborhood graph from rule `neighbors`.
* `.obsp["neighbors"]["connectivities"]` neighborhood graph from rule `neighbors`.
* `.obsp["preprocessing"]` adds 'normalization' and 'log-transformed' info.

## Parameters

You can add these to the config file.

You can add extra arguments to `highly_variable_genes` and `neighbors` as
as `args` the config file.

* normalize.params.raw_counts: H5AD file to raw counts.
* highly_variable_genes.params:
  - batch: batch column to account when getting the HVGs.
  - lineage: if you want to calculate lineage specific genes alone
    or combined with batch.
  - args: more arguments to pass to `scanpy.pp.highly_variable_genes`.
* neighbors.params.args: argumetns for `scanpy.pp.neighbors`.
* pca.params.scale: wether or not to scale counts before calculating
the PCA embedding.

## What happens in each step

**Normalize:** transforms the data using normalize_total and then log1p. It also
makes sure it is a sparse matrix.

**Highly variable genes selection:** highly variable genes are calculated after
using `sc.pp.filter_genes(min_cells=1)` and the `.var` in the unfiltered
object is updated.

**PCA:** if `scale=True` scaling of `.X` will be performed. Then the PCA is
calculated using the highly varaible genes.

**Neighbors:** It will attempt to use the RAPIDS implementation but will default,
if it fails, to the UMAP implementation
[arXiv:1802.03426v3](https://arxiv.org/abs/1802.03426v3).

**UMAP:** It can also use the RAPIDS implementation.

On RAPIDS implementation: whether it is used or not depends on 'os'
in the config file.
