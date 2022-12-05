# Preprocessing

## Input

anndata.Anndata with:

* `.X` containing raw counts
* `.obs` containing all the columns that are defined by `load_data` workflow
* `.var` containing all the columns that are defined by `load_data` workflow

## Output

anndata.Anndata with:

* `.X` normalised and log1 transformed counts from rule `normalize`
* `.var["highly_variable"]` highly variable gene status from rule `highly_variable_genes`
* `.obsm["X_pca"]` PCA from rule `pca`
* `.obsp["distances"]` neighborhood graph from rule `neighbors`
* `.obsp["connectivities"]`neighborhood graph from rule `neighbors`

TODO: save only preprocessing output and add to anndata
