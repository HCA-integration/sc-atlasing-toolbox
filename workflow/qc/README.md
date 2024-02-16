# QC
Quality control of single cell datasets with [sctk autoqc](https://teichlab.github.io/sctk/notebooks/automatic_qc.html) thresholds.
This module provides workflows that quantify the quality of single cell datasets and provide automated thresholds for filtering out low quality cells.
QC metrics and thresholds are visualised in different plots and tables.

## Input

### AnnData file
AnnData file in h5ad or zarr format with the following:

+ raw count matrix under adata.X
+ any categorical or continuous metadata under adata.obs


### Configuration

```yaml
out_dir: test/out
images: test/images

DATASETS:
  dataset_name:
    input:
      qc:
        test: test/input/pbmc68k.h5ad
        test2: test/input/pbmc68k.h5ad
    qc:
      hue:  # colors for colouring joint plots and stratify removed cellls plots
        - phase
        - bulk_labels
      thresholds:  # mapping of file -> custom thresholds to overwrite autoqc thresholds
        test:
          n_counts_max: 500
          n_genes_min: 200
          percent_mito_max: 0.5
```

TODO: explain sctk autoqc threshold naming


## Output

The output of the QC workflow is a set of plots and tables that visualise the quality of the dataset and the thresholds used for filtering cells. The output plots are saved under the images directory from the configuration file.

Outputs include:

* joint scatter plots of QC metrics stratified by `hues` and joint density plots in log scale and regular scale
* barplots of removed cells, stratified by `hues`
* table of thresholds used for filtering cells at the dataset and file level