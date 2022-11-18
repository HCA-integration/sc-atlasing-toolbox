# Subset Data

Subsetting cells for feasible benchmarking

Approaches:

+ Subset cells within samples
+ Subset full samples/datasets
+ others

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
+ `.uns['subset']` subsetting strategy
