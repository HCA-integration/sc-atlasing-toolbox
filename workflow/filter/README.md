## Filtering

This workflow filters the input data based on a set of conditions.
The conditions can be defined in the config under the `filter` module.


```yaml
DATASETS:
  Lee2020:
    input:
      filter: test/input/load_data/harmonize_metadata/Lee2020.zarr
    filter:
      subset: true ## whether to subset the file according to the filters, is False by default
      remove_by_column:
        sample: #column name
        # entries from column to exlucde (will be treated as String)
          - Schulte-Schrepping_C2P01H_d0
          - Schulte-Schrepping_C2P05F_d0
          - Schulte-Schrepping_C2P07H_d0
          - Schulte-Schrepping_C2P10H_d0
          - Schulte-Schrepping_C2P13F_d0
          - Schulte-Schrepping_C2P15H_d0
          - Schulte-Schrepping_C2P16H_d0
          - Schulte-Schrepping_C2P19H_d0
        donor: # column name
        # entries from column to exlucde (will be treated as String)
          - C19-CB-0008 
        disease:
          - influenza
  test:
    input:
      filter: test/input/pbmc68k.h5ad
    filter:
      subset: false
      remove_by_column:
        phase:
          - G1
        is_cd14_mono:
          - true
```

The different conditions in `remove_by_column` are combined with an `AND` operation.
