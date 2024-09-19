# Label transfer

```yaml
DATASETS:
  test:
    input:
      label_transfer:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    label_transfer:
      majority_reference:
        reference_key: bulk_labels
        query_key: louvain
      majority_consensus:
        columns:
          - bulk_labels
          - na_column
          - na_str_column
```

* `majority_reference`: Majority voting to assign reference labels to query clusters.
  * `reference_key`: The key in the reference dataset that contains the labels to be transferred.
  * `query_key`: The key in the query dataset that contains the clusters that will be assigned.
* `majority_consensus`: Majority voting per cell across different clustering assignments.
  * `columns`: The keys in the query dataset containing the same label set.