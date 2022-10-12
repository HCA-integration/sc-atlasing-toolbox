# Integration Methods

## Testing
Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

## Input
AnnData file h5ad with the following:

+ `.uns['dataset']` name of task/dataset
+ `.X` normalised and log-transformed counts (same as `.layers['normcounts']`)
+ `.layers['counts']` raw counts
+ `.layers['normcounts']` normalised and log-transformed counts
+ `.uns['preprocessing']` preprocessing information
  + `.uns['preprocessing']['hvg']`: number of HVGs, 0 means all genes included
  + `.uns['preprocessing']['scaled']`: boolean, whether data is scaled

## Output

Anndata with integrated and unintegrated information:

### Metadata
+ `.uns['dataset']` name of task/dataset
+ `.uns['methods']` methods applied to this object
+ `.uns['integration']` integration specific entries
    + `.uns['integration']['method']` intgration method name
    + `.uns['integration']['label_key']` label used for integration
    + `.uns['integration']['batch_key']` batch used for integration
    + `.uns['integration']['output_type']` output type of method (one of 'knn', 'embed' or 'full')
 

### Unintegrated data
Some slots are dedicated to the unintegrated data and must be available in the integrated anndata object.

+ `.layers['counts']`
+ `.layers['normcounts']`
+ `.raw.obsm['X_pca']`
+ `.uns['preprocessing']` preprocessing parameters as provided

### Integrated output
Different integration methods provide different types of output.
Depending on the output type, the object must contain the following slots:

1. corrected graph (`knn`)
   + `.obsp['connectivities']`, `.obsp['distances']` integrated kNN graph returned by integration method
2. corrected embedding (`embed`)
   + `.obsm['X_emb']` integrated embedding returned by integration method
3. corrected features (`full`)
   + `.X` corrected feature counts returned by integration method
   + `.obsm['X_pca']` PCA on corrected feature counts
