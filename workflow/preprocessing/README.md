# Preprocessing

This module provides rules to processes datasets according to the following steps:

1. normalize counts
2. calculate highly variable genes (and optionally subset to those genes)
3. PCA
4. calculate k-nearest neighbors
5. UMAP dimensionality reduction

Additionally, this module contains rules for plotting a custom embedding (e.g. PCA) or the UMAP.

The rules are defined under `rules/rules.smk` and imported into an end-to-end pipeline in `rules/assemble.smk`

## General input

Each rule can take either a h5ad file or a zarr file as input and will always write the output as a zarr file for more efficient storage.

The input AnnData object must include:

* raw counts (untransformed) in `.X` (by default) or in `.layers` under the key specified with `raw_counts` in the config file (see below for an example).
* `.obs` batch and lineage columns if user wishes to run `highly_variable_genes` with them.

The AnnData files produced by the `load_data` workflow  should work out of the box.

While each step saves only the parts of the AnnData that have changed, the final `assemble` rule will collect the slots that are defined in the config into a single AnnData zarr file.
The assembled file is saved under `config["output_dir"] + "/dataset_name/preprocessed.zarr"`, where `dataset_name` will be replaced by the name that you give your dataset.

In the following example you see all the possible preprocessing slots that the assembled object should contain.

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: adata.h5ad
    preprocessing:
      raw_counts: X  # default
      assemble:
        - counts  # raw counts that are provided as input in .layers['counts']
        - normalize # normalised counts
        - highly_variable_genes
        - pca
        - neighbors
        - umap
```

## Preprocessing steps

Each preprocessing step can be configured with additional parameters in the config, however, if these parameters are not defined, the pipeline will run with default parameters.

### Normalize

Transforms the data using normalize_total and then log-transforms them with `scanpy.pp.log1p`.

**Output**

- `.X` normalised and log1 transformed counts, stored as sparse matrix.
- `.uns["preprocessing"]` containing metadata on the normalisation approach with the following keys:
  - `'normalization'`: normalisation strategy (`'default'` only for now)
  - `'log-transformed'`: whether the counts are log-normalised (always `True` for now)

Note, that any counts stored in `.layers` will be removed and counts in .X will be overwritten.

**Parameters**

- `raw_counts`: Key in .layers or .X itself to specify which matrix gets transformed

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: 'adata.h5ad'
    preprocessing:
      raw_counts: X
      assemble:
        - normalize
```

### Highly variable genes selection

Highly variable genes are calculated after using `scanpy.pp.filter_genes(min_cells=1)` and the `.var` in the unfiltered object is updated.
Note, that the highly variable gene selection will be run on the normalisation output.

**Output**

- `.var["highly_variable"]` boolean column with highly variable gene status.
- `.uns["preprocessing"]["highly_variable_genes"]` containing all the arguments passed to the `scanpy.pp.highly_variable_genes` function (`batch_key` included)

**Parameters**

- `batch`: batch column to account when getting the HVGs.
- `lineage`: if you want to calculate lineage specific genes alone
  or combined with batch.
- `highly_variable_genes`: any additional arguments to pass to the `scanpy.pp.highly_variable_genes` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      highly_variable_genes:
        n_top_genes: 2000
        subset: true  # will subset the feature space of the object, resulting object will be smaller
      assemble:  # only include highly_variable_genes output here, normalized count matrix won't be saved
        - highly_variable_genes
```

### PCA

If `scale=True` scaling of `.X` will be performed, otherwise the counts in `.X` are used directly.
Subsequently, the PCA is calculated using the highly variable genes.
Note, that this step requires the normalized and log-transformed counts in `.X` as well as the highly variable genes information in `.var`.
As part of the pipeline, the `normalization` and `highly_variable_genes` rules will be used as input.

**Output**

- `.obsm["X_pca"]` PCA embedding
- `.uns["preprocessing"]["pca"]` containing all the arguments passed to the `scanpy.pp.pca` function (`use_highly_variable=True` included)
- `.uns["preprocessing"]["scale"]`: whether the matrix was scaled or not before PCA

**Parameters**

- `scale`: wether or not to scale counts before calculating
the PCA embedding.
- `pca`: any additional arguments to pass to the `scanpy.pp.pca` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      pca:
        n_comps: 50
      assemble:
        - pca
```


### K-nearest neighbor graph

It will attempt to use the RAPIDS[^1] implementation but will default, if it fails, to the UMAP implementation
[arXiv:1802.03426v3](https://arxiv.org/abs/1802.03426v3).

By default, the rule will compute the neighbors based on the PCA distances in `.obsm["X_pca"]`, if available, or on `.X` directly otherwise.
In order to use a different representation, `use_rep` must be specified under the additional parameters.

**Output**

- `.obsp["neighbors"]["distances"]` distance matrix
- `.obsp["neighbors"]["connectivities"]` adjacency matrix
- `.uns['preprocessing']['scale']`: whether the matrix was scaled or not before PCA

**Parameters**

- `neighbors`: any additional arguments for the `scanpy.pp.neighbors` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      neighbors:
        use_rep: X_pca
      assemble:
        - neighbors
```

### UMAP

UMAP dimensionality reduction is calculated from the PCA output. It can also use the RAPIDS[^1] implementation.

**Output**

- `.obsm["X_umap"]` UMAP embedding; it becoms .obsm["X_umap_{key}"] if multiple `neighbors_key` are given in the config file.

**Parameters**

- `umap`: any additional arguments for the `scanpy.tl.umap` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      umap:
        neighbors_key: neighbors
      assemble:
        - umap
```

## Assembled output

The `assemble` rule collects the user specified preprocessing steps and saves them in a single anndata.AnnData object (zarr file) with the following slots:

**Output**

- `.layers["counts"]`: raw counts (from input) if `counts` is present under `assemble`
- `.X` and `.layers["normcounts"]`: normalized counts if `normalize` is present under `assemble`
- `.var[["highly_variable", "means", "dispersions", "dispersions_norm", "highly_variable_nbatches", "highly_variable_intersection"]]`: highly variable gene information if `highly_variable_genes` is present under `assemble`
- `.obsm["X_pca"]`: PCA representation if `pca` is present under `assemble`
- `.uns["neighbors"]`, `.obsp["distances"]`, `.obsp["connectivities"]`: kNN graph if `neighbors` is present under `assemble`
- `.obsm["X_umap"]`: UMAP representation if `umap` is present under `assemble`
- `.uns["preprocessing"]` containing the metadata on how the different preprocessing steps were run. See each preprocessing step for more detailed descriptions of the preprocessing metadata.

**Parameters**

- `assemble`: list of preprocessing outputs (including raw counts) to assemble in the final output

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      assemble:
        - counts  # raw counts that are provided as input in .layers['counts']
        - normalize # normalised counts
        - highly_variable_genes
        - pca
        - neighbors
        - umap
```

---

[^1] On RAPIDS implementation: whether it is used or not depends on 'os'
in the config file.
