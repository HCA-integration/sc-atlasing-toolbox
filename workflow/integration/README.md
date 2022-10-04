# Integration Methods

## Testing
Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```commandline
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh --use-conda -c2  # actual run with max 2 cores
```

## Input
AnnData file h5ad with the following:

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

### Unintegrated output
+ `.layers['counts']`
+ `.layers['normcounts']`

### Integrated output
Different integration methods provide different types of output.
These are one of:

1. corrected graph (`knn`)
2. corrected embedding (`embed`)
3. corrected features (`full`)

#### Graph output
For all integration methods.

+`.obsp['connectivities']`: Integrated graph connectivities
+`.obsp['distances']`: Integrated graph distances


#### Embedding output
For most integration methods

+ `.obsm['X_emb']`

#### Feature output
Only for specific integration methods

+ `.layers['corrected_counts']`