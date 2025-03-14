# Label transfer

## Input
Anndata object with:

* `.X`
* label key in `.obs`, specified in `config[DATASETS][<dataset>][label]`
* file name specified in `config[DATASETS][<dataset>][adata_file]`
* label transfer setup specified in `config[DATASETS][<dataset>][label_transfer]`
  * keys: methods e.g. `celltypist, scarches`
  * values: method specific model names e.g. `[Healthy_COVID19_PBMC, Immune_All_Low]`

TODO: Preprocessed anndata

## Output

* TSV file with predictions matched to barcodes
* File of Anndata that has been predicted on


## Methods

* CellTypist
* scArches


## Steps

1. get/prepare model
2. predict from model
