# Metrics API

## Input

### Unintegrated output
+ `.layers['counts']` raw counts
+ `.layers['normcounts']` normalised and log-transformed counts

### Integrated output

One or a list of output types in `.uns['integration']['output_type']` and representations in the corresponding `Anndata`
slots.

1. corrected graph (`knn`)
   + `.obsp['connectivities']` integrated graph connectivities
   + `.obsp['distances']` integrated graph distances
2. corrected embedding (`embed`)
   + `.obsm['X_emb']` integrated embedding
3. corrected features (`full`)
   + `.layers['corrected_counts']` integrated counts


## Output

TSV file containing the metrics per method in long format
e.g.

```tsv
metric method output_type score
nmi    scvi   embed       0.876
asw    scvi   embed       0.876
nmi    scanvi embed       0.605
asw    scanvi embed       0.605
```

Metric names match what is given in the config file.