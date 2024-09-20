# Collect files

Collect slots from files to a single file.
Some modules produce multiple files from a single input file. This module collects these files into a single file.

## Parameters

```yaml
DATASETS:
  SchulteSchrepping2020:
    input:
      collect:
        file_1: test/input/load_data/download/SchulteSchrepping2020.h5ad
        file_2: test/input/load_data/harmonize_metadata/SchulteSchrepping2020_local.zarr
    collect:
      same_slots:
        - X
        - var
      merge_slots:
        - obs
        - obsm
      obs_index_col:
        file_2: barcode
```

* `same_slots`: list of slots that are assumed to be the same in all files. These slots are copied from the first file.
* `merge_slots`: list of slots that are merged from all files. For DataFrame slots in the object, the columns that contain different values will be saved under a renamed column with the file name as suffix. For dictionary slots ()
* `sep`:
* `obs_index_col`: str or mapping of obs column to be used as index for obs. If str, use same column for all files.