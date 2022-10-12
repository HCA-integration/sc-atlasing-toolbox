# Metrics API

## Testing
Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n
bash test/run_test.sh -c
```

## Input

### Unintegrated data
+ `.layers['counts']` raw counts
+ `.layers['normcounts']` normalised and log-transformed counts
+ `.raw` integrated anndata object with all processed slots

### Integrated data

One or a list of output types in `.uns['integration']['output_type']` and representations in the corresponding `Anndata`
slots.

1. corrected graph (`knn`)
   + `.obsp['connectivities']`, `.obsp['distances']` integrated kNN graph returned by integration method
2. corrected embedding (`embed`)
   + `.obsm['X_emb']` integrated embedding returned by integration method
3. corrected features (`full`)
   + `.X` corrected feature counts returned by integration method
   + `.obsm['X_pca']` PCA on corrected feature counts


## Output

TSV file containing the metrics per method in long format with the following columns

+ `metric`: name of metric as specified in params.tsv
+ `method`: integration method
+ `output_type`: one of `knn`, `embed` or `full`, the output type that an integration method returns
+ `metric_type`: one of `bio_conservation` or `batch_correction`
+ `score`: value of metric. NA if metric is not or cannot be computed

e.g.
```tsv
metric      method   output_type   metric_type        score
nmi         scvi     embed         bio_conservation   0.876
asw_label   scvi     embed         bio_conservation   0.876
nmi         scanvi   embed         bio_conservation   0.605
asw_label   scanvi   embed         bio_conservation   0.605
```