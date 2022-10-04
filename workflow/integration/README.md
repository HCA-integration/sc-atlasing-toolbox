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

+ .X
+ .layers['counts']
+ .layers['normcounts']
+ .uns['preprocessing'] preprocessing information
  + hvg: number of HVGs
  + scaled: boolean

## Output

Anndata with integrated and unintegrated information:

### Metadata
+ .uns['dataset'] name of task/dataset
+ .uns['methods'] methods applied to this object
+ .uns['integration'] integration specific entries
    + .uns['integration']['method'] intgration method name
    + .uns['integration']['label_key'] label used for integration
    + .uns['integration']['batch_key'] batch used for integration
    + .uns['integration']['output_type'] output type of method (one of 'knn', 'embed' or 'full')

### Unintegrated output
+ .layers['counts']
+ .layers['normcounts']

### Integrated output

#### Graph output
For all integration methods.

+`.obsp['connectivities']`: Integrated graph connectivities
+`.obsp['distances']`: Integrated graph distances


#### Embedding output
For most integration methods

+ `adata.obsm['X_emb']`

#### Feature output
Only for specific integration methods

+ `adata.layers['integrated_counts']`