# Unimodal Batch Integration

## Input

### AnnData file
AnnData file in h5ad or zarr format with the following:

+ for methods using counts: count matrices (both must be present)
  + raw counts under either `.X`, `.raw.X` or `.layers`
  + normalised and log-transformed counts under either `.X`, `.raw.X` or `.layers`
+ for methods using PCA embedding:
  + PCA embedding under `.obsm['X_pca']`
+ for methods kNN graph:
  + graph connectivities under `.obsp['connectivities']`
  + graph distances under `.obsp['distances']`
  + kNN graph metadata under `.uns['neighbors']` (optional, will speed up unintegrated)

### Configuration

```yaml
DATASETS:
  dataset_name:
    input:
      integration: adata.h5ad
    integration:
      raw_counts: layers/counts # where to find the raw counts in AnnData hiearchy
      norm_counts: layers/normcounts # where to find the normalised counts in AnnData hiearchy
      methods: # list methods to be run
         bbknn:
         combat:
         scanorama:
         unintegrated:
         scgen:
          n_epochs: 10  # hyperparamters for the method
```

## Output

### Metadata
+ `.uns['dataset']` name of task/dataset
+ `.uns['methods']` methods applied to this object
+ `.uns['integration']` integration specific entries
    + `.uns['integration']['method']` intgration method name
    + `.uns['integration']['label_key']` label used for integration
    + `.uns['integration']['batch_key']` batch used for integration
    + `.uns['integration']['output_type']` output type of method (one of 'knn', 'embed' or 'full')


### Integrated output
The direct method output is under `{out_dir}/integration/dataset{dataset}/file_id~{file_id}/batch~{batch}/method~{method}--hyperparams~{hyperparams}--label~{label}/adata.zarr`
Different integration methods provide different types of output.
Depending on the output type, the object must contain the following slots:

1. corrected graph (`knn`)
   + `.obsp['connectivities']`, `.obsp['distances']` integrated kNN graph returned by integration method
2. corrected embedding (`embed`)
   + `.obsm['X_emb']` integrated embedding returned by integration method
3. corrected features (`full`)
   + `.X` corrected feature counts returned by integration method
   + `.obsm['X_pca']` PCA on corrected feature counts

### Processed integrated output
The integrated output files require further processing to be used in downstream analysis.
Files with computed embeddings (for full feature outputs) and kNN graphs (for full feature and embedding outputs) are stored under `{out_dir}/integration/dataset{dataset}/file_id~{file_id}/batch~{batch}/method~{method}--hyperparams~{hyperparams}--label~{label}--output_type~{output_type}.zarr`.

## Testing
Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

The script will call the pipeline and run a test script.
If the input files don't yet exist, the test script might throw a `FileNotFoundError`.
You can ignore it, if you haven't yet executed the pipeline completely and therefore don't yet have all outputs.

### Test model loading
Some integration methods store pytorch model.
Below are tests if these models can be loaded.

```
conda run --live-stream -n scarches python test/run_scarches_model_loading.py 
conda run --live-stream -n scvi-tools python test/run_scvi-tools_model_loading.py 
```
This requires the `scarches` and `scvi-tools` conda environments to be installed, which are defined under `envs/`. 