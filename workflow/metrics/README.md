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
   + `.obsp['connectivities']` integrated graph connectivities from embedding
   + `.obsp['distances']` integrated graph distances from embedding
3. corrected features (`full`)
   + `.layers['corrected_counts']` integrated counts
   + `.obsm['X_emb']` integrated embedding from corrected counts
   + `.obsp['connectivities']` integrated graph connectivities from embedding
   + `.obsp['distances']` integrated graph distances from embedding


## Output

TSV file containing the metrics per method in long format with the following columns

+ metric: name of metric as specified in params.tsv
+ method: integration method
+ output_type: one of `knn`, `embed` or `full`, the output type that an integration method returns
+ metric_type: one of `bio_conservation` or `batch_correction`
+ score: value of metric. NA if metric is not or cannot be computed

e.g.
```tsv
metric      method   output_type   metric_type        score
nmi         scvi     embed         bio_conservation   0.876
asw_label   scvi     embed         bio_conservation   0.876
nmi         scanvi   embed         bio_conservation   0.605
asw_label   scanvi   embed         bio_conservation   0.605
```